// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// NavigationLayer.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_DETECTOR_NAVIGATIONLAYER_H
#define ACTS_DETECTOR_NAVIGATIONLAYER_H

class MsgStream;

// Geometry module
#include "ACTS/Layers/Layer.hpp"
#include "ACTS/Utilities/BinnedArray.hpp"
#include "ACTS/Utilities/Definitions.hpp"

namespace Acts {
    
  class Surface;
  class BinUtility;
    
  /**
   @class NavigationLayer
  
   Class to be used for gaps in Volumes as a navigational link.
   Navigation Layers have a surface representation, but should usually never be
   propagated to.
      
   
   */

 class NavigationLayer : public Layer {
        
      public:
        /** Factory Constructor - the surface representation is given by pointer (ownership passed)
              - spacer layer if needed  */
        static LayerPtr create(Surface* sRepresentation, double thickness=0.)
        { return LayerPtr(new NavigationLayer(sRepresentation, thickness)); }

        /** Clone with a shift - only cloning that is allowed */
        LayerPtr cloneWithShift(const Transform3D& shift) const override;

        /** Destructor*/
        virtual ~NavigationLayer();
        
        /** The binning position method - as default the center is given, but may be overloaded */
        virtual Vector3D binningPosition(BinningValue bValue) const override;
        
        /** Copy Constructor - fobidden */
        NavigationLayer(const NavigationLayer&) = delete;
                                              
        /** Assignment operator - forbidden */
        NavigationLayer& operator=(const NavigationLayer&) = delete;
                    
        /** Transforms the layer into a Surface representation for extrapolation */
        const Surface& surfaceRepresentation() const override;
        
        /** isOnLayer() method, using isOnSurface() with Layer specific tolerance */
        bool isOnLayer(const Vector3D& gp, const BoundaryCheck& bcheck = BoundaryCheck(true)) const override;
        
        /** Boolean check method if layer has material:
           - checks if any of the layer surfaces has material:
           - can be approach surfaces or layer surface */
        bool hasMaterial() const override;

        /** Boolean check method if layer has sensitive surfaces */
        bool hasSensitive() const override;
        
    protected:
        /** Default Constructor*/
        NavigationLayer(){}
      
        /** Constructor - the surface representation is given by pointer (ownership passed)
            - spacer layer if needed  */
        NavigationLayer(Surface* surfaceRepresentation, 
                        double thickness);    
    
        Surface*  m_surfaceRepresentation;       //!< for the navigation Volume the surface is a private member */
      
  };

  inline const Surface&  NavigationLayer::surfaceRepresentation() const 
  { 
    return (*m_surfaceRepresentation); 
  }  
  
  inline Vector3D NavigationLayer::binningPosition(BinningValue bValue) const
  {
      return m_surfaceRepresentation->binningPosition(bValue);
  }  
  
  inline bool NavigationLayer::hasMaterial() const
  {
      return false;
  }

  inline bool NavigationLayer::hasSensitive() const
  {
      return false;
  }
  
} // end of namespace

#endif // ACTS_DETECTOR_NAVIGATIONLAYER_H
