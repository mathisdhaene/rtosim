/* -------------------------------------------------------------------------- *
 * Copyright (c) 2010-2016 C. Pizzolato, M. Reggiani                          *
 *                                                                            *
 * Licensed under the Apache License, Version 2.0 (the "License");            *
 * you may not use this file except in compliance with the License.           *
 * You may obtain a copy of the License at:                                   *
 * http://www.apache.org/licenses/LICENSE-2.0                                 *
 *                                                                            *
 * Unless required by applicable law or agreed to in writing, software        *
 * distributed under the License is distributed on an "AS IS" BASIS,          *
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.   *
 * See the License for the specific language governing permissions and        *
 * limitations under the License.                                             *
 * -------------------------------------------------------------------------- */

#include "rtosim/rtosim.h"
using namespace rtosim;

#include <iostream>
using std::cout;
using std::endl;
#include <string>
using std::string;
#include <Simbody.h>

void printAuthors() {
    cout << "Real-time OpenSim inverse kinematics" << endl;
    cout << "Authors: Claudio Pizzolato <c.pizzolato@griffith.edu.au>" << endl;
    cout << "         Monica Reggiani <monica.reggiani@unipd.it>" << endl << endl;
}

void printHelp() {
    printAuthors();
    auto progName(SimTK::Pathname::getThisExecutablePath());
    bool isAbsolute;
    string dir, filename, ext;
    SimTK::Pathname::deconstructPathname(progName, isAbsolute, dir, filename, ext);

    cout << "Option              Argument         Action / Notes\n";
    cout << "------              --------         --------------\n";
    cout << "-h                                   Print the command-line options for " << filename << ".\n";
    cout << "--model             ModelFilename    Specify the name of the osim model file for the investigation.\n";
    cout << "--trc               TrcFilename      Specify the name of the trc file to be used.\n";
    cout << "--task-set          TaskSetFilename  Specify the name of the XML TaskSet file containing the marker weights to be used.\n";
    cout << "--fc                CutoffFrequency  Specify the name of lowpass cutoff frequency to filter IK data.\n";
    cout << "-j                  IK threads       Specify the number of IK threads to be used.\n";
    cout << "-a                  Accuracy         Specify the IK solver accuracy.\n";
    cout << "--output            OutputDir        Specify the output directory.\n";
    cout << "--push-frequency    PushFrequency    Specify the frequency to which the trajectories are read from the storage file.\n";
    cout << "-v                                   Show visualiser.\n";
}

int main(int argc, char* argv[]) {

    // DEBUG: Print arguments
    std::cout << "DEBUG: Argument Count: " << argc << std::endl;
    for (int i = 0; i < argc; ++i) {
        std::cout << "DEBUG: Arg[" << i << "]: " << argv[i] << std::endl;
    }

    ProgramOptionsParser po(argc, argv);
    if (po.exists("-h") || po.empty()) {
        printHelp();
        exit(EXIT_SUCCESS);
    }

    string osimModelFilename;
    if (po.exists("--model"))
        osimModelFilename = po.getParameter("--model");
    else {
        printHelp();
        exit(EXIT_SUCCESS);
    }
    std::cout << "DEBUG: Loaded model file: " << osimModelFilename << std::endl;

    string trcTrialFilename;
    if (po.exists("--trc"))
        trcTrialFilename = po.getParameter("--trc");
    else {
        printHelp();
        exit(EXIT_SUCCESS);
    }
    std::cout << "DEBUG: Loaded TRC file: " << trcTrialFilename << std::endl;

    string ikTaskFilename;
    if (po.exists("--task-set"))
        ikTaskFilename = po.getParameter("--task-set");
    else {
        printHelp();
        exit(EXIT_SUCCESS);
    }
    std::cout << "DEBUG: Loaded IK Task Set: " << ikTaskFilename << std::endl;

    double fc(8);
    if (po.exists("--fc"))
        fc = po.getParameter<double>("--fc");
    std::cout << "DEBUG: Filter Cutoff Frequency: " << fc << std::endl;

    unsigned nThreads(1);
    if (po.exists("-j"))
        nThreads = po.getParameter<unsigned>("-j");
    std::cout << "DEBUG: Number of IK Threads: " << nThreads << std::endl;

    double solverAccuracy(1e-5);
    if (po.exists("-a"))
        solverAccuracy = po.getParameter<double>("-a");
    std::cout << "DEBUG: Solver Accuracy: " << solverAccuracy << std::endl;

    string resultDir("Output");
    if (po.exists("--output"))
        resultDir = po.getParameter("--output");
    std::cout << "DEBUG: Output Directory: " << resultDir << std::endl;

    double pushFrequency(-1);
    if (po.exists("--push-frequency"))
        pushFrequency = po.getParameter<double>("--push-frequency");
    std::cout << "DEBUG: Push Frequency: " << pushFrequency << std::endl;

    bool showVisualiser(false);
    if (po.exists("-v"))
        showVisualiser = true;
    std::cout << "DEBUG: Visualiser Enabled: " << showVisualiser << std::endl;

    resultDir = rtosim::FileSystem::getAbsolutePath(resultDir);
    rtosim::FileSystem::createDirectory(resultDir);
    string stopWatchResultDir(resultDir);

    // DEBUG: Queues and Barriers
    std::cout << "DEBUG: Setting up queues and synchronization barriers..." << std::endl;

    rtosim::MarkerSetQueue markerSetQueue;
    rtosim::GeneralisedCoordinatesQueue generalisedCoordinatesQueue, filteredGeneralisedCoordinatesQueue;

    rtb::Concurrency::Latch doneWithSubscriptions;
    rtb::Concurrency::Latch doneWithExecution;

    auto coordNames = getCoordinateNamesFromModel(osimModelFilename);
    rtosim::GeneralisedCoordinatesStateSpace gcFilt(fc, coordNames.size());

    // DEBUG: Setting up threads...
    std::cout << "DEBUG: Setting up threads..." << std::endl;

    rtosim::MarkersFromTrc markersFromTrc(markerSetQueue, doneWithSubscriptions, doneWithExecution, osimModelFilename, trcTrialFilename, false);
    markersFromTrc.setOutputFrequency(pushFrequency);

    rtosim::QueueToInverseKinematics inverseKinematics(markerSetQueue, generalisedCoordinatesQueue, doneWithSubscriptions, doneWithExecution, osimModelFilename, nThreads, ikTaskFilename, solverAccuracy);

    rtosim::QueueAdapter<rtosim::GeneralisedCoordinatesQueue, rtosim::GeneralisedCoordinatesQueue, rtosim::GeneralisedCoordinatesStateSpace>
    gcQueueAdaptor(generalisedCoordinatesQueue, filteredGeneralisedCoordinatesQueue, doneWithSubscriptions, doneWithExecution, gcFilt);

    rtosim::QueueToFileLogger<rtosim::GeneralisedCoordinatesData> filteredIkLogger(filteredGeneralisedCoordinatesQueue, doneWithSubscriptions, doneWithExecution, coordNames, resultDir, "filtered_ik_from_file", "sto");

    rtosim::QueueToFileLogger<rtosim::GeneralisedCoordinatesData> rawIkLogger(generalisedCoordinatesQueue, doneWithSubscriptions, doneWithExecution, coordNames, resultDir, "raw_ik_from_file", "sto");

    rtosim::FrameCounter<rtosim::GeneralisedCoordinatesQueue> ikFrameCounter(generalisedCoordinatesQueue, "time-ik-throughput");

    rtosim::TimeDifference<rtosim::GeneralisedCoordinatesQueue, rtosim::GeneralisedCoordinatesQueue>
    gcQueueAdaptorTimeDifference(generalisedCoordinatesQueue, filteredGeneralisedCoordinatesQueue, doneWithSubscriptions, doneWithExecution);

    rtosim::TimeDifference<rtosim::MarkerSetQueue, rtosim::GeneralisedCoordinatesQueue>
    ikTimeDifference(markerSetQueue, generalisedCoordinatesQueue, doneWithSubscriptions, doneWithExecution);

    doneWithSubscriptions.setCount(7);
    doneWithExecution.setCount(7);

    // DEBUG: Launching threads
    std::cout << "DEBUG: About to launch threads, visualiser = " << showVisualiser << std::endl;

    if (showVisualiser) {
        rtosim::StateVisualiser visualiser(generalisedCoordinatesQueue, osimModelFilename);
        rtosim::QueuesSync::launchThreads(markersFromTrc, inverseKinematics, gcQueueAdaptor, filteredIkLogger, gcQueueAdaptorTimeDifference, ikTimeDifference, rawIkLogger, ikFrameCounter, visualiser);
    } else {
        rtosim::QueuesSync::launchThreads(markersFromTrc, inverseKinematics, gcQueueAdaptor, filteredIkLogger, gcQueueAdaptorTimeDifference, ikTimeDifference, rawIkLogger, ikFrameCounter);
    }

    std::cout << "DEBUG: Threads completed, post-processing..." << std::endl;

    auto stopWatches = inverseKinematics.getProcessingTimes();
    rtosim::StopWatch combinedSW("time-ikparallel-processing");
    for (auto& s : stopWatches)
        combinedSW += s;
    combinedSW.print(stopWatchResultDir);
    ikFrameCounter.getProcessingTimes().print(stopWatchResultDir);
    ikTimeDifference.getWallClockDifference().print(stopWatchResultDir + "/time-markerqueue-to-jointangles.txt");
    gcQueueAdaptorTimeDifference.getWallClockDifference().print(stopWatchResultDir + "/time-jointangles-to-filteredjointangles.txt");

    std::cout << "DEBUG: Program completed successfully." << std::endl;
    return 0;
}
