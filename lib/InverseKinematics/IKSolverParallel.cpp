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
#include <OpenSim/Common/Constant.h>
#include <OpenSim/Common/Set.h>
#include <OpenSim/Simulation/MarkersReference.h>
#include <OpenSim/Tools/IKCoordinateTask.h>
#include <SimTKcommon.h>

#include <OpenSim/Simulation/InverseKinematicsSolver.h>
#include <memory>
using std::unique_ptr;
#include <limits>
#include <iostream>
#include <algorithm>

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
        coordinateTaskConfigs_.clear();
        for (size_t i(0); i < static_cast<size_t>(ikTaskSet.getSize()); ++i) {
            const auto& task = ikTaskSet.get(i);
            auto& taskMutable = const_cast<OpenSim::IKTask&>(task);
            if (!task.getApply()) {
                continue;
            }

            std::string currentMarkerName(task.getName());
            auto it = markerWeights_.find(currentMarkerName);
            if (it != markerWeights_.end()) {
                markerWeights_[task.getName()] = taskMutable.getWeight();
            }

            const auto* coordTask = dynamic_cast<const OpenSim::IKCoordinateTask*>(&task);
            if (coordTask) {
                CoordinateTaskConfig cfg;
                cfg.name = task.getName();
                cfg.weight = taskMutable.getWeight();
                cfg.valueType = static_cast<int>(coordTask->getValueType());
                cfg.value = coordTask->getValue();
                coordinateTaskConfigs_.push_back(cfg);
            }
        }

        std::cerr << "[IKSolverParallel] Loaded coordinate tasks: "
                  << coordinateTaskConfigs_.size() << std::endl;
    }

    void IKSolverParallel::setInverseKinematicsTaskSet(const std::string& ikTaskSetFilename) {
        OpenSim::IKTaskSet ikTaskSet(ikTaskSetFilename);
        setInverseKinematicsTaskSet(ikTaskSet);
    }

    void IKSolverParallel::pushState(const SimTK::State& s) {
        GeneralisedCoordinatesData currentData(nCoordinates_);
        std::vector<double> q(nCoordinates_);
        model_.realizePosition(s);
        for (unsigned i(0); i < nCoordinates_; ++i) {
            q[i] = model_.getCoordinateSet().get(i).getValue(s);
        }
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
        coordinateRefs.clear();
        for (const auto& cfg : coordinateTaskConfigs_) {
            if (!model_.getCoordinateSet().contains(cfg.name)) {
                continue;
            }
            double targetValue = cfg.value;
            if (cfg.valueType == static_cast<int>(OpenSim::IKCoordinateTask::DefaultValue)) {
                targetValue = model_.getCoordinateSet().get(cfg.name).getDefaultValue();
            } else if (cfg.valueType == static_cast<int>(OpenSim::IKCoordinateTask::FromFile)) {
                targetValue = model_.getCoordinateSet().get(cfg.name).getValue(s);
            }
            OpenSim::Constant constantValue(targetValue);
            OpenSim::CoordinateReference coordRef(cfg.name, constantValue);
            coordRef.setWeight(cfg.weight);
            coordinateRefs.push_back(coordRef);
        }

        SimTK::State defaultState(s);
        unsigned ct = 0;
        unsigned failCount = 0;

        if (parityMode_) {
            OpenSim::InverseKinematicsSolver ikSolver(model_, *markerReference, coordinateRefs, contraintWeight_);
            ikSolver.setAccuracy(sovlerAccuracy_);
            ikSolver.setAdvanceTimeFromReference(false);
            ikSolver.assemble(s);
            defaultState = s;

            while (localRunCondition) {
                if (!markerReference->isEndOfData()) {
                    double currentTime = inputThreadPoolJobs_.front().time;
                    s.updTime() = currentTime;
                    try {
                        ikSolver.track(s);
                    } catch (const std::exception& e) {
                        ++failCount;
                        std::cerr << "[IKSolverParallel] parity track() failed at frame #" << ct
                                  << " time=" << currentTime
                                  << " : " << e.what() << std::endl;
                        s = defaultState;
                        try { ikSolver.assemble(s); } catch (...) {}
                    } catch (...) {
                        ++failCount;
                        std::cerr << "[IKSolverParallel] parity track() failed at frame #" << ct
                                  << " time=" << currentTime
                                  << " : unknown exception" << std::endl;
                        s = defaultState;
                        try { ikSolver.assemble(s); } catch (...) {}
                    }

                    if ((ct % 50) == 0 && !isWithinRom(s)) {
                        std::cerr << "[IKSolverParallel] Frame #" << ct
                                  << " has coordinates outside model ranges." << std::endl;
                    }
                    if ((ct % 50) == 0) {
                        SimTK::Array_<double> markerErrs;
                        ikSolver.computeCurrentMarkerErrors(markerErrs);
                        double errSqSum = 0.0;
                        double maxErr = 0.0;
                        for (int i = 0; i < markerErrs.size(); ++i) {
                            const double e = markerErrs[i];
                            errSqSum += e * e;
                            maxErr = std::max(maxErr, e);
                        }
                        const int usedMarkers = markerErrs.size();
                        const double rmsErr = usedMarkers > 0 ? std::sqrt(errSqSum / usedMarkers)
                                                               : std::numeric_limits<double>::quiet_NaN();
                        std::cerr << "[IKSolverParallel] Frame #" << ct
                                  << " time=" << currentTime
                                  << " usedMarkers=" << usedMarkers
                                  << " disabledMarkers=0"
                                  << " rmsErr(m)=" << rmsErr
                                  << " maxErr(m)=" << maxErr
                                  << " failures=" << failCount << std::endl;
                    }

                    pushState(s);
                    defaultState = s;
                    ++ct;
                } else {
                    localRunCondition = false;
                    outputGeneralisedCoordinatesQueue_.push(rtosim::EndOfData::get<GeneralisedCoordinatesFrame>());
                }
            }
        } else {
            while (localRunCondition) {
                if (!markerReference->isEndOfData()) {
                OpenSim::Set<OpenSim::MarkerWeight> frameWeights;
                SimTK::Array_<SimTK::Vec3> markerVals;
                markerReference->getValues(s, markerVals);
                // `getValues()` pops the next frame and updates the internal time.
                // Use that time for both the state and the per-frame table below.
                const double currentTime = markerReference->getCurrentTime();
                s.updTime() = currentTime;

                for (int i = 0; i < markerVals.size(); ++i) {
                    // Preserve task-set marker weights and only disable truly invalid data.
                    double weight = markerWeights_.at(markerNames_[i]);
                    if (std::isnan(markerVals[i][0]) || std::isnan(markerVals[i][1]) || std::isnan(markerVals[i][2])) {
                        weight = 0.0;
                    }
                    // Many TRC pipelines encode missing markers as (0,0,0) rather than NaN.
                    // Treat those as occluded to avoid pulling the model to the origin.
                    if (markerVals[i].normSqr() < 1e-18) {
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
                } catch (const std::exception& e) {
                    ++failCount;
                    std::cerr << "[IKSolverParallel] track() failed at frame #" << ct
                              << " time=" << currentTime
                              << " : " << e.what() << std::endl;
                    s = defaultState;
                } catch (...) {
                    ++failCount;
                    std::cerr << "[IKSolverParallel] track() failed at frame #" << ct
                              << " time=" << currentTime
                              << " : unknown exception" << std::endl;
                    s = defaultState;
                }
                if ((ct % 50) == 0 && !isWithinRom(s)) {
                    std::cerr << "[IKSolverParallel] Frame #" << ct
                              << " has coordinates outside model ranges." << std::endl;
                }

                // Debug: marker fit diagnostics for this frame.
                model_.realizePosition(s);
                double errSqSum = 0.0;
                double maxErr = 0.0;
                int usedMarkers = 0;
                int disabledMarkers = 0;
                for (int i = 0; i < markerVals.size(); ++i) {
                    double weight = markerWeights_.at(markerNames_[i]);
                    if (std::isnan(markerVals[i][0]) || std::isnan(markerVals[i][1]) || std::isnan(markerVals[i][2]) ||
                        markerVals[i].normSqr() < 1e-18) {
                        weight = 0.0;
                    }
                    if (weight <= 0.0) {
                        ++disabledMarkers;
                        continue;
                    }

                    const auto& modelMarker = model_.getMarkerSet().get(markerNames_[i]);
                    const SimTK::Vec3 modelPos = modelMarker.getLocationInGround(s);
                    const double err = (modelPos - markerVals[i]).norm();
                    errSqSum += err * err;
                    maxErr = std::max(maxErr, err);
                    ++usedMarkers;
                }
                const double rmsErr = usedMarkers > 0 ? std::sqrt(errSqSum / usedMarkers) : std::numeric_limits<double>::quiet_NaN();
                if ((ct % 50) == 0) {
                    std::cerr << "[IKSolverParallel] Frame #" << ct
                              << " time=" << currentTime
                              << " usedMarkers=" << usedMarkers
                              << " disabledMarkers=" << disabledMarkers
                              << " rmsErr(m)=" << rmsErr
                              << " maxErr(m)=" << maxErr
                              << " failures=" << failCount << std::endl;
                }

                std::vector<double> qVals(nCoordinates_);
                for (unsigned i = 0; i < nCoordinates_; ++i) {
                    qVals[i] = model_.getCoordinateSet().get(i).getValue(s);
                }
                std::cerr << "[IKSolverParallel] Q values frame #" << ct << ": ";
                for (int i = 0; i < static_cast<int>(qVals.size()); ++i) {
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
        }

        doneWithExecution_.wait();
    }

    IKSolverParallel::~IKSolverParallel() {
#ifdef RTOSIM_DEBUG
        cout << " IKSolver " << std::this_thread::get_id() << " is closing" << endl;
#endif
    }
}
