#!/usr/bin/perl

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# 
# Description:
# ----------- 
#    This script is a driver routine used to plot contour plots
#    of big number of datasets.
#                                                radostin simitev
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

#---------------------------------------------------
# variables
#---------------------------------------------------
$fortran         = 'Yokoi-plots-rs';
$F_batch         = 'Yokoi.in';
$mima_txt        = '_mima.txt';
$_idl_vis_lev    = '/home/staff2/dynamo/07_DRS/03_Visualisation/02_ContourPlots/IDL_scripts/o_Contour';
$format          = 2;
$_idl_batch_file = 'idl.batch';
$levelnum        = 21;
$idl             = '/maths/rsi/idl/bin/idl';


#----------------------------------------------------
# 1     Read the input arguments to the perl script:
# - data directory
# - plots destination directory
# - dataset name
# - numbers to plot
#----------------------------------------------------
my $path     = shift @ARGV;
my $dest_dir = shift @ARGV;
my $dataset  = shift @ARGV;
my @numbers  =       @ARGV;
my @params   = ('EMF', 'alphaBz', 'betaJz', 'gammaOz' );

if($path eq '.') {
    $path = `pwd`;
    chomp($path);
}

chdir($path);

#----------------------------------------------------
# Make the input file for the fortran program
#----------------------------------------------------
my $f_batch = '';
foreach $i (@numbers) {
   $run_name=$dataset.".".$i;
   $f_batch .= "$run_name \n"
}

#--------------------------------------------------
# Make Fortran batch. Run Fortran to make a cut
#--------------------------------------------------
#...make fortran batch file
open(FILE,"> $path/$F_batch")
    or die "Cannot write to $path/$F_batch: $!\n";
print FILE $f_batch;
close(FILE);

#...execute the fortran cutting program
system("$fortran");

#...move files to the plots folder
system("mkdir -p $dest_dir");
system("mv idl.dims $dest_dir");
system("mv *.x $dest_dir");
system("mv *.y $dest_dir");
system("mv *.z $dest_dir");
system("mv *_mima.txt $dest_dir");
system("mv Steady* $dest_dir");
system("mv *_rms.txt $dest_dir");

chdir($dest_dir);

#---------------------------------------------------
# Find min & max of data
#---------------------------------------------------
foreach $par (@params) {

   open(FILE,"$par$mima_txt") or die "Cannot open $par$mima_txt: $!\n";
   while($mm=<FILE>) {
       chomp($mm);
       ($min,$max) = split(' ',$mm);
       if($lmax > $max) { $max = $lmax }
       if($lmin < $min) { $min = $lmin }
   }
   close(FILE);

   print " \n" ;
   print "abs. Minimum: $min\n";
   print "abs. Maximum: $max\n";
   print " \n" ;


   #---------------------------------------------------
   #  Create IDL batch file and execute idl < batch
   #---------------------------------------------------
   open(FILE,"idl.dims") or die "Cannot read idl.dims : $!\n";
   chomp(($ncosf,$lmaxf,$nanglef,$m0,$eta) = <FILE>);
   close(FILE);
   @parms =('');

   push(@parms,"abs_min = $min");
   push(@parms,"abs_max = $max");
   push(@parms,"levelnum = $levelnum");
   push(@parms,"ncosf = $ncosf ");
   push(@parms,"lmaxf = $lmaxf");
   push(@parms,"nanglef = $nanglef");
   push(@parms,"m0 = $m0");
   push(@parms,"eta = $eta");

   $r_parm_list = \@parms;

   # create IDL batch files
   my $rsh_file = _make_idl_batch_file($_idl_vis_lev,$par,$format,$r_parm_list);
   open(FILE,"> $_idl_batch_file")
      or die "Cannot write to $_idl_batch_file: $!\n";

   print FILE $rsh_file;
   close(FILE);

   # call IDL
   print "plot of $par ...  "; 
   system("$idl < $_idl_batch_file >& /dev/null");
   # rename plots
   system("mv plot.ps $par.plot.ps");
   print "(done)\n";

}

# ----- private methods -----

sub _make_idl_batch_file {
    my $idl_script = shift;
    my $param = shift;
    my $format = shift;
    my $r_parm_list = shift;

    $a = '';
    $a   .= "tmp_x = \'$param.x\'\n";
    $a   .= "tmp_y = \'$param.y\'\n";
    $a   .= "tmp_z = \'$param.z\'\n";
    $a   .= "fo = $format\n";
    $a   .= "cut_type = 1\n";
    $a   .= "cut_piece = 7\n";
    $a   .= "part = 1\n";
    $a   .= "field = 3\n";
    $a   .= "imean = 0\n";
    $a   .= "phi = 0.0\n";
    foreach (@$r_parm_list) { $a .= "$_\n"; }
    $a   .= ".run $idl_script\n";
    $a   .= "exit\n";

    return($a);
}
# vim: tabstop=3:softtabstop=3:shiftwidth=3:expandtab
