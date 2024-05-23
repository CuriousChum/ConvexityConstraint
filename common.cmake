# --- !! TO USER : Change these two paths !! ---

# On JMM's macbook
SET(Conda_ConvexityConstraint_Env_Dir "/home/vboxuser/miniconda3/envs/ConvexityConstraint/")
SET(HamiltonFastMarching_Source_Dir "/home/vboxuser/programming/cpp/HamiltonFastMarching")

# The above two paths correspond to:
# - A ConvexityConstraint conda environment, containing the cgal, nlopt, and eigen packages, using the provided yaml file. This can be done with the following command : conda env create -f ConvexityConstraint.yaml
# - The HamiltonFastMarching repository contents, which can be been downloaded at : https://github.com/Mirebeau/HamiltonFastMarching

# ------------------ PC version (not working yet) ------------------------

# On JMM's PC
#SET(Conda_ConvexityConstraint_Env_Dir "C:/Users/jmmir/miniconda3/envs/ConvexityConstraint")
#SET(HamiltonFastMarching_Source_Dir "C:/Users/jmmir/Documents/GitHub/HamiltonFastMarching")

#set(CMAKE_PREFIX_PATH ${CondaDir} ${CondaDir}/Library/lib/CGAL)

#Set(CgalLib  "${CondaDir}/Library/lib/CGAL-vc140-mt-4.9.1.lib")   # CGAL
#Set(NloptLib "${CondaDir}/Library/lib/nlopt.lib")  # NLOPT
#Set(GmpLib   "${CondaDir}/Library/lib/gmp.lib")    # GMP


# ------------------- Mac version ----------------------

# General setup
set(CondaDir ${Conda_ConvexityConstraint_Env_Dir})

# --- Set library paths ---
set(CMAKE_PREFIX_PATH ${CondaDir})
find_package(CGAL)

Set(CgalLib  "${CondaDir}/lib/libCGAL.so")   # CGAL
Set(NloptLib "${CondaDir}/lib/libnlopt.so")  # NLOPT
Set(GmpLib   "${CondaDir}/lib/libgmp.so")    # GMP


# --- Set header paths ---
Set(ExternalHeadersDir
	"${CondaDir}/include/eigen3"   # Eigen
#	"/Users/mirebeau/Dropbox/Programmes/2022/eigen-3.2.3"
	"${CondaDir}/include"          # CGAL and NLOPT
	"${HamiltonFastMarching_Source_Dir}/JMM_CPPLibs"  # A few personnal classes and macros used
	)

