#!/bin/bash
# Set defaults
TEST=serial
BUILD_TYPE=RELEASE

# Usage instructions.
function usage()
{
	echo "`basename $0` [-t <name>] [-B <type>]"
	echo "`basename $0` -?"
	echo " <name> The name of the test. Can be one of 'serial'(default), 'parallel', 'serialDD' or 'parallelDD'."
	echo " <type> The build type. can be one of 'DEBUG' or 'RELEASE'(default)"
}

function do_compile()
{
	mkdir -p ${BUILDDIR}
	cd ${BUILDDIR}
	cmake -DCMAKE_BUILD_TYPE=$BUILD_TYPE \
	      -DBUILD_DOCS=OFF \
	      -DBUILD_UTILS=OFF \
	      -DBUILD_TESTS=OFF \
	      -DCOMP=$COMP \
	      -DMPI=$MPI \
	      -DCMAKE_INSTALL_PREFIX=${INSTALLDIR}  ..
	make clean
	make && make install|| exit
	cd ${BASEDIR}
	rm -rf ${BUILDDIR}
}

function do_run_test()
{
	# Go to the test folder and copy all the initial
	# files there
	cd ${TESTDIR}
	cp ${BASEDIR}/examples/states/BENCHMARK/* .
	# Start the run
	(
	./jobscript $TEST

	eb=`tail -n1 e035p1t2r100000m1p5test.2.eb|awk '{printf "%s", $2}'`
	ek=`tail -n1 e035p1t2r100000m1p5test.2.ek|awk '{printf "%s", $2}'`

	echo '==========================================================='
	echo " Actual eb     |      Expected eb "
	echo " $eb |    0.913453561D+04 "
	test "$eb" = "0.913453561D+04" && echo "Magnetic energy test passed!"|| (echo "Magnetic energy test failed!"; exit 1)

	echo '==========================================================='
	echo " Actual ek     |      Expected ek "
	echo " $ek |    0.448474975D+03 "
	test "$ek" = "0.448474975D+03" && echo "Kinetic energy test passed!" || (echo "Kinetic energy test failed!"; exit 1)
	) &
	cd $BASEDIR
}

function do_plot()
{
	cd ${TESTDIR}
	# Create a live plot that we will use to track the state of the 
	# simulation with
	cat - << EOF > plot_energies.gpi
	reset
	set datafile fortran
	plot 'e035p1t2r100000m1p5test.2.ek' u 1:2 axes x1y1 w lp t 'kinetic now', \
	     '../test-BENCHMARK/e035p1t2r100000m1p5test.2.ek' u 1:2 axes x1y1 w l t 'kinetic ref', \
	     'e035p1t2r100000m1p5test.2.eb' u 1:2 axes x1y2 w lp t 'magnetic now', \
	     '../test-BENCHMARK/e035p1t2r100000m1p5test.2.eb' u 1:2 axes x1y2 w l t 'magnetic ref'
	pause 10
	reread
EOF

	# Sleep until the relevant file was created.
	until test -f e035p1t2r100000m1p5test.2.ek
	do
	   sleep 1
	done
	gnuplot plot_energies.gpi 
	rm drs.lock
	cd $BASEDIR
}

# Parse arguments.
while getopts ":t:B" opt
do
   case $opt in
      t)
         TEST=$OPTARG
         ;;
      B)
         BUILD_TYPE=$OPTARG
         ;;
      \?)
         echo "Invalid option: -$OPTARG" >&2
	 usage
         exit 1
         ;;
      :)
         echo "Option -$OPTARG requires an argument." >&2
	 usage
         exit 2
         ;;
   esac
done

BASEDIR=${PWD}
TESTDIR=${BASEDIR}/test-$TEST
BUILDDIR=${BASEDIR}/build-$TEST
INSTALLDIR=${TESTDIR}

# Set compilation flags.
case $TEST in
	'serial')
		MPI=OFF
		COMP=OFF
		;;
	'parallel')
		MPI=ON
		COMP=OFF
		;;
	'serialDD')
		MPI=OFF
		COMP=ON
		;;
	'parallelDD')
		MPI=ON
		COMP=ON
		;;
esac


# Clean up the test folder.
# This needs to be done before compilation and installation or we lose everything.
if [[ -d ${TESTDIR} ]]
then
	echo " WARNING: This run will destroy the contents of ${TESTDIR}!"
	echo -n " ??? Should I continue? (Y/N) "
	read cont
	test "$cont"=="Y" || exit
	rm -rf ${TESTDIR}
fi
mkdir ${TESTDIR}
# Build the binaries and install everything to the run folder.
do_compile

# Run the test. We do it in a subshell so we can then track the state.
do_run_test 

# Make a real-time plot.
do_plot

