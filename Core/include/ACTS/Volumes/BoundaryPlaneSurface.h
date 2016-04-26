///////////////////////////////////////////////////////////////////
// BoundaryPlaneSurface.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_VOLUMES_BOUNDARYPLANESURFACE_H
#define ACTS_VOLUMES_BOUNDARYPLANESURFACE_H 1

// Geometry module
#include "ACTS/Surfaces/PlaneSurface.h"
#include "ACTS/Volumes/BoundarySurface.h"
// EventData module
#include "ACTS/Utilities/PropDirection.h"
// Core module
#include "ACTS/Utilities/AlgebraDefinitions.h"

namespace Acts {

    class Volume;

  /**
   @class BoundaryPlaneSurface

   BoundaryPlaneSurface description inside the tracking realm,
   it extends the PlaneSurface description to make a surface being a boundary of a
   Acts::Volume (used for all volume shapes).
   It inherits from BoundarySurface to get the interface of boundaries.

   @author Andreas.Salzburger@cern.ch
  */

  template <class T> class BoundaryPlaneSurface :
                               virtual public BoundarySurface<T>, public PlaneSurface {

     typedef std::shared_ptr<const T> VolumePtr;
     typedef BinnedArray< VolumePtr > VolumeArray;

    public:
     /** Default Constructor - needed for pool and inherited classes */
     BoundaryPlaneSurface() :
       BoundarySurface<T>(),
       PlaneSurface()
     {}

     /** Copy constructor */
     BoundaryPlaneSurface(const BoundaryPlaneSurface<T>& bps) :
       BoundarySurface<T>(bps),
       PlaneSurface(bps)
     {}

     /** Constructor for a Boundary with exact two Volumes attached to it*/
     BoundaryPlaneSurface(const T* inside, const T* outside, const PlaneSurface& psf) :
       BoundarySurface<T>(inside, outside),
       PlaneSurface(psf)
     {}

     /** Constructor for a Boundary with two VolumeArrays attached to it*/
     BoundaryPlaneSurface(std::shared_ptr<const VolumeArray> insideArray, std::shared_ptr<const VolumeArray> outsideArray, const PlaneSurface& psf) :
       BoundarySurface<T>(insideArray, outsideArray),
       PlaneSurface(psf)
     {}

     /** Copy constructor with a shift */
     BoundaryPlaneSurface(const T* inside, const T* outside, const PlaneSurface& psf, const Transform3D& tr) :
       BoundarySurface<T>(inside,outside),
       PlaneSurface(psf,tr)
     {}

     /** The Surface Representation of this */
     const Surface& surfaceRepresentation() const override;

     /**Virtual Destructor*/
     virtual ~BoundaryPlaneSurface(){}

     /**Assignment operator - addignment of boundary surfaces are forbidden*/
     BoundaryPlaneSurface& operator=(const BoundaryPlaneSurface& vol);

  };

template <class T> inline const Surface& BoundaryPlaneSurface<T>::surfaceRepresentation() const { return *this; }

} // end of namespace Acts

#endif // ACTS_VOLUMES_BOUNDARYPLANESURFACE_H
