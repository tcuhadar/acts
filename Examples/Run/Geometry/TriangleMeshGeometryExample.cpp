// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <Acts/Geometry/GeometryContext.hpp>
#include <Acts/Geometry/TrackingGeometry.hpp>
#include <Acts/Visualization/GeometryView.hpp>


#include <memory>
#include <string>
#include <vector>

#include "ACTFW/Detector/IBaseDetector.hpp"
#include "ACTFW/Framework/AlgorithmContext.hpp"
#include "ACTFW/Framework/IContextDecorator.hpp"
#include "ACTFW/Framework/WhiteBoard.hpp"
#include "ACTFW/Geometry/CommonGeometry.hpp"
#include "ACTFW/GenericDetector/GenericDetector.hpp"
#include "ACTFW/Options/CommonOptions.hpp"
#include "ACTFW/Utilities/Options.hpp"
#include "ACTFW/Utilities/Paths.hpp"
#include <ACTFW/Geometry/TriangleMesh.hpp>

#include "ACTFW/Plugins/Obj/ObjTrackingGeometryWriter.hpp"
#include "ACTFW/Plugins/Obj/ObjWriterOptions.hpp"

int processGeometry(int argc, char* argv[], FW::IBaseDetector& detector) {
  // setup and parse options
  auto desc = FW::Options::makeDefaultOptions();
  FW::Options::addSequencerOptions(desc);
  FW::Options::addGeometryOptions(desc);
  FW::Options::addMaterialOptions(desc);
  FW::Options::addObjWriterOptions(desc);
  FW::Options::addOutputOptions(desc);

  // Add specific options for this geometry
  detector.addOptions(desc);
  auto vm = FW::Options::parse(desc, argc, argv);
  if (vm.empty()) {
    return EXIT_FAILURE;
  }

  // Now read the standard options
  auto logLevel = FW::Options::readLogLevel(vm);
  auto nEvents = FW::Options::readSequencerConfig(vm).events;

  // The geometry, material and decoration
  auto geometry = FW::Geometry::build(vm, detector);
  auto tGeometry = geometry.first;
  auto contextDecorators = geometry.second;

  // The detectors
  read_strings subDetectors = vm["geo-detector-volume"].as<read_strings>();

  auto volumeLogLevel =
      Acts::Logging::Level(vm["geo-volume-loglevel"].as<size_t>());

  for (size_t ievt = 0; ievt < nEvents; ++ievt) {
    // Setup the event and algorithm context
    FW::WhiteBoard eventStore(
        Acts::getDefaultLogger("EventStore#" + std::to_string(ievt), logLevel));
    size_t ialg = 0;

    // The geometry context
    FW::AlgorithmContext context(ialg, ievt, eventStore);

    /// Decorate the context
    for (auto& cdr : contextDecorators) {
      if (cdr->decorate(context) != FW::ProcessCode::SUCCESS)
        throw std::runtime_error("Failed to decorate event context");
    }

    std::string geoContextStr = "";
    if (contextDecorators.size() > 0) {
      // We need indeed a context object
      if (nEvents > 1) {
        geoContextStr = "_geoContext" + std::to_string(ievt);
      }
    }

    // Convert to TriangleMesh
    auto world = tGeometry->highestTrackingVolume();
    if (world) {
      auto cfg = FW::Options::readObjTrackingGeometryWriterConfig(
          vm, "ObjTrackingGeometryWriter", volumeLogLevel);
      cfg.containerView.triangulate = true;
      cfg.volumeView.triangulate = true;
      cfg.sensitiveView.triangulate = true;
      cfg.passiveView.triangulate = true;
      cfg.gridView.triangulate = true;

      Acts::TriangleMesh triangleMesh(cfg.outputPrecision, cfg.outputScalor);

      Acts::GeometryView::drawTrackingVolume(
        triangleMesh, *world, context.geoContext, cfg.containerView,
        cfg.volumeView, cfg.sensitiveView, cfg.passiveView, cfg.gridView);

      std::cout << triangleMesh << std::endl;
    }



    // // ---------------------------------------------------------------------------------
    // // Output directory
    // std::string outputDir = vm["output-dir"].template as<std::string>();

    // // OBJ output
    // if (vm["output-obj"].as<bool>()) {
    //   // Configure the tracking geometry writer
    //   auto tgObjWriterConfig = FW::Options::readObjTrackingGeometryWriterConfig(
    //       vm, "ObjTrackingGeometryWriter", volumeLogLevel);
    //   auto tgObjWriter = std::make_shared<FW::Obj::ObjTrackingGeometryWriter>(
    //       tgObjWriterConfig);
    //   // Write the tracking geometry object
    //   tgObjWriter->write(context, *tGeometry);
    // }
  }

  return 0;
}

/// @brief main executable
///
/// @param argc The argument count
/// @param argv The argument list
int main(int argc, char* argv[]) {
  // --------------------------------------------------------------------------------
  GenericDetector detector;
  // now process it
  return processGeometry(argc, argv, detector);
}
