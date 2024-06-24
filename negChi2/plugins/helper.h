#ifndef PhysicsTools_BParkingNano_helpers
#define PhysicsTools_BParkingNano_helpers

#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/GeometryCommonDetAlgo/interface/Measurement1D.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/GeometryCommonDetAlgo/interface/GlobalError.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "DataFormats/GeometryVector/interface/PV3DBase.h"
#include "Math/LorentzVector.h"

// for the fit
#include "RecoVertex/KinematicFit/interface/TwoTrackMassKinematicConstraint.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicTree.h"
#include "RecoVertex/KinematicFit/interface/KinematicConstrainedVertexFitter.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicConstraint.h"
#include "RecoVertex/KinematicFitPrimitives/interface/MultiTrackKinematicConstraint.h"

// 4 vectors
#include "TLorentzVector.h"
#include "TVector3.h" // for boost vector

// for the cov matrix correction
#include "TMatrixDSym.h"
#include "TVectorD.h"
#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"

// for the reco function
#include <tuple>

#include <vector>
#include <algorithm>
#include <limits>
#include <memory>

typedef std::vector<reco::TransientTrack> TransientTrackCollection;


///////////////////////////////////////////////////////////////////////////////////
// Fix the track covariance matrix to be pos. def.


inline reco::Track correctCovMat(const reco::Track *tk, double delta){

// Parameters associated to the 5D curvilinear covariance matrix: 
// (qoverp, lambda, phi, dxy, dsz) 
// Defined as:
//   qoverp = q / abs(p) = signed inverse of momentum [1/GeV] 
//   lambda = pi/2 - polar angle at the given point 
//   phi = azimuth angle at the given point 
//   dxy = -vx*sin(phi) + vy*cos(phi) [cm] 
//   dsz = vz*cos(lambda) - (vx*cos(phi)+vy*sin(phi))*sin(lambda) [cm] 

    unsigned int i;
    unsigned int j;
    double min_eig = 1;

    // Get the original covariance matrix
    reco::TrackBase::CovarianceMatrix cov = tk->covariance();

    // Define a TMatrixDSym of the same shape as the old cov matrix.
    // Sym -> you only have to give one dimension, it will be a symm matrix
    // of shape (cov.kRows, cov.kRows)
    TMatrixDSym newCov(cov.kRows);

    // loop over old cov matrix 
    for (i = 0; i < cov.kRows; i++) {
        for (j = 0; j < cov.kRows; j++) {
            // change nan or inf values to 1e-6 
            if (std::isnan(cov(i,j)) || std::isinf(cov(i,j)))
                cov(i,j) = 1e-6;
            // fill new covariacne matrix
            newCov(i,j) = cov(i,j);
        }
    }

    // Define a vector of size cov.kRows
    TVectorD eig(cov.kRows);
    // Fill it with the egienvalues of the newCov
    newCov.EigenVectors(eig);

    // loop over eigenvalues and find the minimal eigenvalue :)
    for (i = 0; i < cov.kRows; i++)
        if (eig(i) < min_eig)
            min_eig = eig(i);

    // If the minimum eigenvalue is less than zero, then subtract it from the diagonal and add `delta`.
    if (min_eig < 0) {
        for (i = 0; i < cov.kRows; i++)
            cov(i,i) -= min_eig - delta;
    }

    return reco::Track(tk->chi2(), tk->ndof(), tk->referencePoint(), tk->momentum(), tk->charge(), cov, tk->algo(), (reco::TrackBase::TrackQuality) tk->qualityMask());
    }

///////////////////////////////////////////////////////////////////////////////////
// Fix the track covariance matrix to be pos. def.
inline reco::Track fixTrack(const reco::TrackRef& tk)
{
    reco::Track t = reco::Track(*tk);
    return correctCovMat(&t, 1e-8);
}


///////////////////////////////////////////////////////////////////////////////////
// Vertex Fit 
inline RefCountedKinematicTree vertexFit(std::vector<RefCountedKinematicParticle> toFit, ParticleMass massConstr, bool applyConstr)
{

  //define fitter
  KinematicConstrainedVertexFitter fitter;

  //define constraint
  MultiTrackKinematicConstraint* constr = new TwoTrackMassKinematicConstraint(massConstr);

  RefCountedKinematicTree fitTree;

  if(applyConstr) {
  // perform the fit
  fitTree = fitter.fit(toFit, constr);
  }
  else{
  fitTree = fitter.fit(toFit);
  }

  return fitTree;

}
#endif
