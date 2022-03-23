# --- !! TO USER : Change these two paths !! ---

SET(Conda_ConvexityConstraint_Env_Dir "/Users/mirebeau/opt/miniconda3/envs/ConvexityConstraint/")
SET(HamiltonFastMarching_Source_Dir "/Users/mirebeau/Dropbox/Programmes/Github/HamiltonFastMarching")

# The above two paths correspond to:
# - A ConvexityConstraint conda environment, containing the cgal, nlopt, and eigen packages, using the provided yaml file. This can be done with the following command : conda env create -f ConvexityConstraint.yaml
# - The HamiltonFastMarching repository contents, which can be been downloaded at : https://github.com/Mirebeau/HamiltonFastMarching

# ------------------------------------------


# General setup
set(CondaDir ${Conda_ConvexityConstraint_Env_Dir})

# --- Set library paths ---
set(CMAKE_PREFIX_PATH ${CondaDir})
find_package(CGAL)

Set(CgalLib  "${CondaDir}/lib/libCGAL.dylib")   # CGAL
Set(NloptLib "${CondaDir}/lib/libnlopt.dylib")  # NLOPT
Set(GmpLib   "${CondaDir}/lib/libgmp.dylib")    # GMP


# --- Set header paths ---
Set(ExternalHeadersDir
	"${CondaDir}/include/eigen3"   # Eigen
	"${CondaDir}/include"          # CGAL and NLOPT
	"${HamiltonFastMarching_Source_Dir}/JMM_CPPLibs"  # A few personnal classes and macros used
	)

