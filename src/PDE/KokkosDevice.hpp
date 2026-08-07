#ifndef KokkosDevice_h
#define KokkosDevice_h

#include "Types.hpp"
#include "Exception.hpp"
#include "Fields.hpp"
#include "UnsMesh.hpp"
#include "Integrate/Quadrature.hpp"
#include "Kokkos_Core.hpp"

// Fields storage U(e*nprop+c) is in row major
// Throw error if not detected
#ifndef FIELD_DATA_LAYOUT_AS_FIELD_MAJOR
  #error "FIELD_DATA_LAYOUT must be FIELD_MAJOR. Reconfigure with -DFIELD_DATA_LAYOUT=field_major to fix, or port kernels to EqCompUnk."
#endif

// Max. no. of quadrature pts (11 for NGvol(10))
constexpr std::size_t NQUAD_MAX=14;
// Max. no. of DOFs (10 for DGP2)
constexpr std::size_t NDOF_MAX=10;
// Max. no. of materials supported by device kernels
constexpr std::size_t NMAT_MAX=4;
// Max. no. of conserved components is 3*nmat + 3 + 9*nsld
constexpr std::size_t NCOMP_MAX=50;
// Max. size of combined consv+prim state vector at quad pt
constexpr std::size_t NSTATE_MAX=50;
// Max. size of THINC vol-fraction work array is nmat*rdof
constexpr std::size_t NALSOL_MAX = NMAT_MAX*NDOF_MAX;

using execution_space = Kokkos::DefaultExecutionSpace;
using memory_space = Kokkos::DefaultExecutionSpace::memory_space;

namespace tk {
	// Verify config fits device kernel scratch capacity
	// Throws error if any limit is exceeded
	// Because the thread-private scratch will be corrupted
	// and there is no diagnostic to catch it, so need to exit
	inline void
	checkKokkosCaps( std::size_t nmat,
	                 std::size_t ndof,
	                 std::size_t rdof,
	                 std::size_t ncomp,
	                 std::size_t nprim )
	{
		ErrChk( nmat <= NMAT_MAX, "Max. " + std::to_string(NMAT_MAX) + " materials, got " + std::to_string(nmat)+".");
		ErrChk( ndof <= NDOF_MAX && rdof <= NDOF_MAX, "Max. " + std::to_string(NDOF_MAX) + " DOFs, got {ndof,rdof}={" + std::to_string(ndof)+","+std::to_string(rdof)+"}.");
		ErrChk( ncomp <= NCOMP_MAX, "Max. " + std::to_string(NCOMP_MAX) + " components, got " + std::to_string(ncomp)+".");
		ErrChk( ncomp+nprim <= NSTATE_MAX, "Max. " + std::to_string(NSTATE_MAX) + " state entries, got " + std::to_string(ncomp*nprim)+".");
		ErrChk( nmat*rdof <= NALSOL_MAX, "Max. " + std::to_string(NALSOL_MAX) + " THINC entries, got " + std::to_string(nmat*rdof)+".");
		ErrChk( NGvol(ndof) <= NQUAD_MAX, "Max. " + std::to_string(NQUAD_MAX) + " quadrature points supported.");
	}


	// Persistent device-resident buffers for const-P DG kernels
	// Bundles all device views that would otherwise be allocated by kernels
	// Any instance must be owned by something destroyed before Kokkos::finalize() and released on chare migration
	// Device pointers held are not PUP'd!
	struct KokkosDeviceViews {
		// Time-invariant
		Kokkos::View< std::size_t*, memory_space > solidx;
		Kokkos::View< std::size_t*, memory_space > inpoel;
		Kokkos::View< real*, memory_space > cx;
		Kokkos::View< real*, memory_space > cy;
		Kokkos::View< real*, memory_space > cz;
		Kokkos::View< real*, memory_space > geoElem;

		// Per-call
		Kokkos::View< real*, memory_space > U;
		Kokkos::View< real*, memory_space > P;
		Kokkos::View< real*, memory_space > R;
		Kokkos::View< real*, memory_space > riemannDeriv;

		// Host-provenance of mesh data resident on device
		const std::size_t* src_inpoel = nullptr;
		const real* src_cx = nullptr;
		const real* src_cy = nullptr;
		const real* src_cz = nullptr;
		const real* src_geoElem = nullptr;
		std::size_t src_nelem = 0;
		std::size_t src_npoin = 0;
		std::size_t src_nmat = 0;
		
		// The below logic helps to capture the AMR and signal to reupload the data to the device
		// Mesh generation source (resident-copied state)
		std::size_t src_gen = static_cast< std::size_t > ( -1 );
		// Current mesh generation bumped by owner when mesh changes
		std::size_t gen = 0;
		// Mark resident data as stale without freeing
		void invalidate() noexcept
		{ src_gen = static_cast< std::size_t > ( -1 ); ++gen; }
		// Free every device alloc and reset all residency bookkeeping
		// Must be called from owner pup on migration
		void release() { *this = KokkosDeviceViews(); }
	}; //kdv

        // This is a test to check for the resident mesh if it belongs to the current function call	
	// True if device already has this partition's mesh data at current generation state -> skip the H2D copy
	// TODO:  Get the input parameter from DG indicating if AMR and ALE has happened or not
	inline bool
	meshResident ( KokkosDeviceViews& d,
	               const std::vector< std::size_t >& inpoel,
	               const UnsMesh::Coords& coord,
	               const Fields& geoElem,
	               std::size_t nelem,
	               std::size_t nmat )
	{
		const bool hit = d.src_gen == d.gen
		              && d.src_inpoel == inpoel.data()
		              && d.src_cx == coord[0].data()
		              && d.src_cy == coord[1].data()
		              && d.src_cz == coord[2].data()
		              && d.src_geoElem == geoElem.getPointer()
		              && d.src_nelem == nelem
		              && d.src_npoin == coord[0].size()
		              && d.src_nmat == nmat;
		if (!hit) {
			d.src_gen = d.gen;
			d.src_inpoel = inpoel.data();
			d.src_cx = coord[0].data();
			d.src_cy = coord[1].data();
			d.src_cz = coord[2].data();
			d.src_geoElem = geoElem.getPointer();
			d.src_nelem = nelem;
			d.src_npoin = coord[0].size();
			d.src_nmat = nmat;
		}
		return hit;
	} //mr

	// Ensures persistent device view has requested extent
	// True if reallocated -> contents are uninitialized, needs reupload
	template< typename T >
	bool
	ensureDeviceCapacity( Kokkos::View< T*, memory_space >& view,
	                      const std::string& label,
	                      std::size_t n )
	{
		if(view.extent(0) != n){
			view = Kokkos::View< T*, memory_space >(
				Kokkos::view_alloc( label, Kokkos::WithoutInitializing ), n );
			return true;
		}
		return false;
	} //edc
} //tk

#endif // KokkosDevice_h
