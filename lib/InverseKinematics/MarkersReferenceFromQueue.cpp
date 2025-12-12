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

#include "rtosim/MarkersReferenceFromQueue.h"
#include "rtosim/EndOfData.h"
#include <iostream>
using std::cout;
using std::endl;
#include <vector>
using std::vector;

namespace rtosim{

    MarkersReferenceFromQueue::MarkersReferenceFromQueue(
        ThreadPoolJobs<MarkerSetFrame>& inputMarkerSetFrameQueue,
        const std::vector<std::string>& markerNames,
        const std::vector<double>& markerWeigths) :
        inputMarkerSetFrameQueue_(inputMarkerSetFrameQueue),
        markerWeigths_(markerWeigths),
        time_(std::numeric_limits<double>::infinity()),
        lastQueueTime_(-std::numeric_limits<double>::infinity()){

        for (auto& n : markerNames)
            markerNames_.push_back(n);
    }

    bool MarkersReferenceFromQueue::isEndOfData() const {

        MarkerSetFrame frame(inputMarkerSetFrameQueue_.front());
        return rtosim::EndOfData::isEod(frame);
    }

	void MarkersReferenceFromQueue::purgeCurrentFrame() {
		MarkerSetFrame frame = inputMarkerSetFrameQueue_.pop();  // ✅ Retire le frame courant
		lastQueueTime_ = time_;  // Met à jour le dernier temps utilisé

		// ⚠️ On suppose que la queue a encore des frames. Sinon, time_ devient infini.
		try {
		    time_ = inputMarkerSetFrameQueue_.front().time;  // On regarde juste le prochain time
		}
		catch (...) {
		    // Si jamais on ne peut pas accéder à front(), on met time_ à infini pour indiquer la fin.
		    time_ = std::numeric_limits<double>::infinity();  
		}

		std::cerr << "[MarkersReferenceFromQueue] Purged frame, next time_: " << time_ << std::endl;
	}


	void MarkersReferenceFromQueue::getValues(const SimTK::State &s, SimTK::Array_<SimTK::Vec3> &values) const {
		values.clear();
		MarkerSetFrame frame;
		double currentTime(inputMarkerSetFrameQueue_.front().time);
		if (currentTime < lastQueueTime_)
		    frame = getPastFrame(currentTime);
		else
		    frame = getFrameFromQueue();

		if (rtosim::EndOfData::isEod(frame))
		    for (auto& m : frame.data)
		        m.setOccluded(true);
		for (auto& marker : frame.data)
		    values.push_back(marker.getCoordinates());

		time_ = frame.time;

		std::cerr << "[MarkersReferenceFromQueue] Updated time_: " << time_ << std::endl;
	}



    MarkerSetFrame MarkersReferenceFromQueue::getPastFrame(double time) const {

        auto it(std::lower_bound(pastFrames_.begin(), pastFrames_.end(), time, [](
            MarkerSetFrame& lhs,
            double t){ return lhs.time < t; }));
        return *it;
    }

    MarkerSetFrame MarkersReferenceFromQueue::getFrameFromQueue() const {

        MarkerSetFrame frame(inputMarkerSetFrameQueue_.pop());

        lastQueueTime_ = frame.time;
        pastFrames_.push_back(frame);
        if (pastFrames_.size() < MaxFramesToStore)
            pastFrames_.pop_front();
        return frame;
    }

    const SimTK::Array_<std::string>& MarkersReferenceFromQueue::getNames() const {

        return markerNames_;
    }

    void MarkersReferenceFromQueue::getWeights(const SimTK::State &s, SimTK::Array_<double> &weights) const {

        weights = markerWeigths_;
    }
	void MarkersReferenceFromQueue::setMarkerIsOccluded(int index, bool occluded) {
		if (!pastFrames_.empty()) {
		    MarkerSetFrame& lastFrame = pastFrames_.back();  // Récupère le frame en cours
		    if (index >= 0 && index < static_cast<int>(lastFrame.data.size())) {
		        lastFrame.data[index].setOccluded(occluded);
		    }
    }
}

}
