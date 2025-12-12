/* -------------------------------------------------------------------------- *
 * Copyright (c) 2010-2016 C. Pizzolato, M. Reggiani                          *
 * Licensed under the Apache License, Version 2.0                             *
 * -------------------------------------------------------------------------- */

#include "rtosim/IKSolverParallel.h"
#include "rtosim/MarkersReferenceFromQueue.h"
#include "rtosim/ArrayConverter.h"
#include "rtosim/EndOfData.h"
#include "rtosim/GeneralisedCoordinatesData.h"
#include "rtosim/queue/GeneralisedCoordinatesQueue.h"
using rtosim::GeneralisedCoordinatesData;
using rtosim::GeneralisedCoordinatesFrame;

//#include <OpenSim/OpenSim.h>
#include <OpenSim/Simulation/Model/Model.h>
#include <OpenSim/Simulation/InverseKinematicsSolver.h>
#include <OpenSim/Common/TimeSeriesTable.h>
#include <OpenSim/Common/Set.h>
#include <OpenSim/Simulation/MarkersReference.h>
#include <SimTKcommon.h>

#include <OpenSim/Simulation/InverseKinematicsSolver.h>
#include <memory>
using std::unique_ptr;
#include <limits>
#include <iostream>

namespace rtosim{

    IKSolverParallel::IKSolverParallel(
        ThreadPoolJobs<MarkerSetFrame>& inputThreadPoolJobs,
        IKoutputs<rtosim::GeneralisedCoordinatesFrame>& outputGeneralisedCoordinatesQueue,
        rtb::Concurrency::Latch& doneWithSubscriptions,
        rtb::Concurrency::Latch& doneWithExecution,
        const std::string& osimModelFilename,
        double solverAccuracy,
        double contraintWeight
    ) :
        inputThreadPoolJobs_(inputThreadPoolJobs),
        outputGeneralisedCoordinatesQueue_(outputGeneralisedCoordinatesQueue),
        doneWithSubscriptions_(doneWithSubscriptions),
        doneWithExecution_(doneWithExecution),
        osimModelFilename_(osimModelFilename),
        model_(osimModelFilename),
        sovlerAccuracy_(solverAccuracy),
        contraintWeight_(contraintWeight) {

        OpenSim::Array<std::string> markerNamesArray, coordinateNamesArray;
        const_cast<OpenSim::MarkerSet&>(model_.getMarkerSet()).getMarkerNames(markerNamesArray);
        rtosim::ArrayConverter::toStdVector(markerNamesArray, markerNames_);
        nMarkers_ = markerNames_.size();

        model_.getCoordinateSet().getNames(coordinateNamesArray);
        rtosim::ArrayConverter::toStdVector(coordinateNamesArray, coordinateNames_);
        nCoordinates_ = model_.getNumCoordinates();

        for (auto it : markerNames_)
            markerWeights_.insert(std::make_pair(it, 1)); //init weights to 1
    }

    void IKSolverParallel::setInverseKinematicsTaskSet(const OpenSim::IKTaskSet& ikTaskSet) {
        for (size_t i(0); i < static_cast<size_t>(ikTaskSet.getSize()); ++i) {
            std::string currentMarkerName(ikTaskSet.get(i).getName());
            auto it = markerWeights_.find(currentMarkerName);
            if (it != markerWeights_.end() && ikTaskSet.get(i).getApply()) {
                markerWeights_[ikTaskSet.get(i).getName()] = ikTaskSet.get(i).getWeight();
            }
        }
    }

    void IKSolverParallel::pushState(const SimTK::State& s) {
        GeneralisedCoordinatesData currentData(nCoordinates_);
        std::vector<double> q(nCoordinates_);
        SimTK::Vector stateQ(s.getQ());
        for (unsigned i(0); i < nCoordinates_; ++i)
            q[i] = stateQ[i];
        currentData.setQ(q);
        outputGeneralisedCoordinatesQueue_.push({ s.getTime(), currentData });
    }

    bool IKSolverParallel::isWithinRom(const SimTK::State& s) const {
        bool isInRom(true);
        auto q(s.getQ());
        for (unsigned i(0); i < nCoordinates_; ++i) {
            auto rangeMax(model_.getCoordinateSet().get(i).getRangeMax());
            auto rangeMin(model_.getCoordinateSet().get(i).getRangeMin());
            if (q[i] > rangeMax || q[i] < rangeMin) {
                isInRom = false;
                std::cerr << coordinateNames_[i] << " is outside its range of motion" << std::endl;
            }
        }
        return isInRom;
    }

    void IKSolverParallel::operator()() {
        SimTK::State s = model_.initSystem();
        bool localRunCondition(true);
        std::vector<double> sortedMarkerWeights;
        for (auto it : markerNames_)
            sortedMarkerWeights.push_back(markerWeights_[it]);

        std::cerr << "[IKSolverParallel] Number of markers: " << nMarkers_ << std::endl;

        unique_ptr<MarkersReferenceFromQueue> markerReference(new MarkersReferenceFromQueue(inputThreadPoolJobs_, markerNames_, sortedMarkerWeights));

        OpenSim::Set<OpenSim::MarkerWeight> osimMarkerWeights;
        for (auto it : markerNames_) {
            osimMarkerWeights.adoptAndAppend(new OpenSim::MarkerWeight(it, markerWeights_[it]));
        }
        markerReference->setMarkerWeightSet(osimMarkerWeights);

        doneWithSubscriptions_.wait();
        SimTK::Array_<OpenSim::CoordinateReference> coordinateRefs;

        // 🔵 Initial Solveur pour Assemblage
        OpenSim::InverseKinematicsSolver ikSolver(model_, *markerReference, coordinateRefs, contraintWeight_);
        ikSolver.setAccuracy(sovlerAccuracy_);
        ikSolver.assemble(s);

        SimTK::State defaultState(s);
        pushState(s);

        unsigned ct = 0;
        while (localRunCondition) {
            if (!markerReference->isEndOfData()) {
                double currentTime = markerReference->getCurrentTime();
                s.updTime() = currentTime;

                OpenSim::Set<OpenSim::MarkerWeight> frameWeights;
                SimTK::Array_<SimTK::Vec3> markerVals;
                markerReference->getValues(s, markerVals);

                for (int i = 0; i < markerVals.size(); ++i) {
                    double weight = 1.0;
                    if (std::isnan(markerVals[i][0]) || std::isnan(markerVals[i][1]) || std::isnan(markerVals[i][2])) {
                        weight = 0.0;
                    }
                    frameWeights.adoptAndAppend(new OpenSim::MarkerWeight(markerNames_[i], weight));
                }

                OpenSim::TimeSeriesTable_<SimTK::Vec3> markerTable;
                SimTK::RowVector_<SimTK::Vec3> markerRow(static_cast<int>(markerVals.size()));
                for (int i = 0; i < markerVals.size(); ++i) {
                    markerRow[i] = markerVals[i];
                }
                markerTable.appendRow(currentTime, markerRow);
                markerTable.setColumnLabels(markerNames_);

                unique_ptr<OpenSim::MarkersReference> dynamicMarkersRef(
                    new OpenSim::MarkersReference(markerTable, frameWeights, OpenSim::Units::Meters));

                OpenSim::InverseKinematicsSolver ikSolverTemp(model_, *dynamicMarkersRef, coordinateRefs, contraintWeight_);
                ikSolverTemp.setAccuracy(sovlerAccuracy_);
                ikSolverTemp.assemble(s);

                try {
                    ikSolverTemp.track(s);
                } catch (...) {
                    s = defaultState;
                }

                SimTK::Vector qVals = s.getQ();
                std::cerr << "[IKSolverParallel] Q values frame #" << ct << ": ";
                for (int i = 0; i < qVals.size(); ++i) {
                    std::cerr << qVals[i] << " ";
                }
                std::cerr << std::endl;

                pushState(s);
                defaultState = s;
                ++ct;

                // FIXED: Removed purgeCurrentFrame() to prevent double pop
                // markerReference->purgeCurrentFrame();

            } else {
                localRunCondition = false;
                outputGeneralisedCoordinatesQueue_.push(rtosim::EndOfData::get<GeneralisedCoordinatesFrame>());
            }
        }

        doneWithExecution_.wait();
    }

    IKSolverParallel::~IKSolverParallel() {
#ifdef RTOSIM_DEBUG
        cout << " IKSolver " << std::this_thread::get_id() << " is closing" << endl;
#endif
    }
}

