#!/usr/bin/perl -w
use strict;
use Cwd;
use Getopt::Long;

my $base_base_dir = cwd();

my $base_options_file = "$base_base_dir/Moocho.opt";
# Can be overridden by --extra-args!

# mpi_go = "mpirun -machinefile ~/bin/linux/machinelist-sahp6556 -np"
my $mpi_go = "";

my $exe = "$base_base_dir/NLPThyraEpetraAdvDiffReactOpt.exe";

my $len_x = 1.0;
my $len_y = 1.0;
my $local_nx = 1;
my $local_ny = -1;
my $global_ny = -1;

my $reaction_rate = 1.0;
my $beta = 1.0;

my $starting_num_procs = 1;
my $max_num_procs = 1;

my $num_procs_per_cluster=-1;

my $starting_model_np = 1;
my $max_model_np = 1;
my $incr_model_np = 1;

my $use_amesos = 0;
my $amesos_params_file = "";
my $use_aztecoo = 0;
my $aztecoo_params_file = "";
my $use_belos = 0;
my $belos_params_file = "";

my $starting_block_size = 1;
my $max_block_size = 1;
my $incr_block_size = 1;

my $solve_fwd_prob = 0;
my $x0 = 0.0001;
my $p0 = 1.0;
my $p_solu_scale = 1.001;
my $fwd_init_options_file = "";
my $fwd_options_file = "";

my $do_analytic = 0;
my $inv_options_file = "";

my $do_finite_diff = 0;
my $fd_step_len = -1.0;
my $fd_options_file = "";

my $use_finite_diff_obj = 0;
my $use_finite_diff_con = 0;

my $do_black_box = 0;
my $bb_fd_step_len = 1e-8;
my $bb_options_file = "";

my $use_direct = 0;
my $fwd_newton_tol = 1e-10;

my $use_default_inv_obj = 0;

my $extra_args = "";

GetOptions(
  "mpi-go=s"                => \$mpi_go,
  "exe=s"                   => \$exe,
  "len-x=f"                 => \$len_x,
  "len-y=f"                 => \$len_y,
  "local-nx=i"              => \$local_nx,
  "local-ny=i"              => \$local_ny,
  "global-ny=i"             => \$global_ny,
  "reaction-rate=f"         => \$reaction_rate,
  "beta=f"                  => \$beta,
  "starting-num-procs=i"    => \$starting_num_procs,
  "max-num-procs=i"         => \$max_num_procs,
  "num-procs-per-cluster=s" => \$num_procs_per_cluster,
  "starting-model-np=i"     => \$starting_model_np,
  "max-model-np=i"          => \$max_model_np,
  "incr-model-np=i"         => \$incr_model_np,
  "use-amesos!"             => \$use_amesos,
  "amesos-params-file=s"    => \$amesos_params_file,
  "use-aztecoo!"            => \$use_aztecoo,
  "aztecoo-params-file=s"   => \$aztecoo_params_file,
  "use-belos!"              => \$use_belos,
  "belos-params-file=s"     => \$belos_params_file,
  "starting-block-size=i"   => \$starting_block_size,
  "max-block-size=i"        => \$max_block_size,
  "incr-block-size=i"       => \$incr_block_size,
  "base-options-file=s"     => \$base_options_file,
  "solve-fwd-prob!"         => \$solve_fwd_prob,
  "x0=f"                    => \$x0,
  "p0=f"                    => \$p0,
  "p-solu-scale=f"          => \$p_solu_scale,
  "fwd-init-options-file=s" => \$fwd_init_options_file,
  "fwd-options-file=s"      => \$fwd_options_file,
  "do-analytic!"            => \$do_analytic,
  "inv-options-file=s"      => \$inv_options_file,
  "do-finite-diff!"         => \$do_finite_diff,
  "fd-step-len=f"           => \$fd_step_len,
  "fd-options-file=s"       => \$fd_options_file,
  "use-finite-diff-obj"     => \$use_finite_diff_obj,
  "use-finite-diff-con"     => \$use_finite_diff_con,
  "do-black-box!"           => \$do_black_box,
  "bb-fd-step-len=f"        => \$bb_fd_step_len,
  "bb-options-file=s"       => \$bb_options_file,
  "use-direct!"             => \$use_direct,
  "fwd-newton-tol=f"        => \$fwd_newton_tol,
  "use-default-inv-obj!"    => \$use_default_inv_obj,
  "extra-args=s"            => \$extra_args
  );

my $do_serial =  ($mpi_go eq "" );

$max_num_procs = $starting_num_procs if($starting_num_procs > $max_num_procs);
$max_model_np = $starting_model_np if($starting_model_np > $max_model_np);
$max_block_size = $starting_block_size if($starting_block_size > $max_block_size);

make_abs_path(\$exe);
make_abs_path(\$base_options_file);
make_abs_path(\$amesos_params_file);
make_abs_path(\$aztecoo_params_file);
make_abs_path(\$belos_params_file);
make_abs_path(\$fwd_init_options_file);
make_abs_path(\$fwd_options_file);
make_abs_path(\$inv_options_file);
make_abs_path(\$fd_options_file);
make_abs_path(\$bb_options_file);

mkchdir("runs");

my $rtn = 0;

my $model_np = $starting_model_np;
my $vary_model_np = ( $model_np < $max_model_np );
for( ; $model_np <= $max_model_np; $model_np += $incr_model_np ) {

  mkchdir("model_np-$model_np") if($vary_model_np);

  my $num_procs = $starting_num_procs;
  my $vary_num_procs = ( $num_procs < $max_num_procs && !$do_serial );
  for( ; $num_procs <= $max_num_procs; $num_procs *= 2 ) {

    mkchdir("num_procs-$num_procs") if($vary_num_procs);

    my $num_procs_per_cluster_used = -1;
    if($num_procs_per_cluster > 0) {
      $num_procs_per_cluster_used = $num_procs_per_cluster;
    }
    else {
      $num_procs_per_cluster_used = $num_procs;
    }

    my $len_y_used = $len_y;

    if($global_ny > 0) {
      $local_ny = $global_ny / $num_procs_per_cluster_used;
    }
    else {
      $len_y_used = $len_y * $num_procs_per_cluster_used;
    }

    my $cmnd_first =
      ( $do_serial ? "" : "$mpi_go $num_procs" )
      ." $exe --moocho-options-file=$base_options_file"
      ." --len-x=$len_x --len-y=$len_y_used --local-nx=$local_nx --local-ny=$local_ny"
      ." --reaction-rate=$reaction_rate --beta=$beta"
      ." --np=$model_np --x0=$x0 --p0=$p0"
      ;

    if($num_procs_per_cluster > 0) {
      $cmnd_first .= " --num-procs-per-cluster=$num_procs_per_cluster";
    }

    my $cmnd_last =
      " $extra_args"
      ;

    if( $use_amesos && $num_procs==1 ) {

      mkchdir("amesos");

      my $cmnd =
        $cmnd_first
        ." ".extra_linear_solver_cmndline("Amesos")
        ." --linear-solver-params-file=\"$amesos_params_file\""
        .$cmnd_last
        ;

      $rtn = 1 if (0!=run_case($cmnd));

      chdir "..";

    }

    if($use_aztecoo) {

      mkchdir("aztecoo");

      my $cmnd =
        $cmnd_first
        ." ".extra_linear_solver_cmndline("AztecOO")
        ." --linear-solver-params-file=\"$aztecoo_params_file\""
        .$cmnd_last
        ;

      $rtn = 1 if (0!=run_case($cmnd));

      chdir "..";

    }

    if($use_belos) {

      mkchdir("belos");

      my $cmnd =
        $cmnd_first
        ." ".extra_linear_solver_cmndline("Belos")
        ." --linear-solver-params-file=\"$belos_params_file\""
        .$cmnd_last
        ;

      $rtn = 1 if (0!=run_case($cmnd,1,0));

      my $block_size = min($starting_block_size,$model_np+1);
      my $loop_max_block_size = min($max_block_size,$model_np+1);
      my $vary_block_size = ( $block_size < $loop_max_block_size );
      for( ; $block_size <= $loop_max_block_size; $block_size += $incr_block_size ) {

        mkchdir("block_size-$block_size") if($vary_block_size);

        my $cmnd =
          $cmnd_first
          ." --lowsf=belos --lowsf-params-file=\"$belos_params_file\""
          ." --lowsf-extra-params=\"<ParameterList><Parameter name=\\\"Block Size\\\" type=\\\"int\\\" value=\\\"$block_size\\\"/></ParameterList>\""
          .$cmnd_last
          ;

        $rtn = 1 if (0!=run_case($cmnd,0,1,($vary_block_size?"../../fwd":"../fwd")));

        chdir("..") if($vary_block_size);

      }

      chdir "..";

    }

    chdir("..") if($vary_num_procs);

  }

  chdir("..") if($vary_model_np);

}

chdir "..";

print "\nFINAL END RESULT: TEST PASSED\n" if $rtn==0;
print "\nFINAL END RESULT: TEST FAILED\n" if $rtn!=0;

exit $rtn;

#
# Subroutines
#

sub min {
  my $a = shift;
  my $b = shift;
  return ( $a < $b ? $a : $b );
}

sub make_abs_path {
  my $file_name_ref = shift;
  return if($$file_name_ref eq "");
  if(substr($$file_name_ref,0,1) eq ".") {
    $$file_name_ref = $base_base_dir . "/" . $$file_name_ref;
  }
}

sub mkchdir {
  my $dir_name = shift;
  mkdir $dir_name, 0777;
  chdir $dir_name;
}

sub extra_linear_solver_cmndline {
  my $solver_name = shift;
  return "--extra-linear-solver-params=\""
    ."<ParameterList>"
      ."<Parameter name=\\\"Linear Solver Type\\\" type=\\\"string\\\" value=\\\"${solver_name}\\\"/>"
    ."</ParameterList>\""
}

sub moocho_thyra_solver_cmndline {
  #
  my $cmnd_hashref = shift;
  #
  my $do_sim = exists $cmnd_hashref->{"do-sim"};
  my $use_direct = $use_direct;
  $use_direct = 0 if(exists $cmnd_hashref->{"use-first-order"});
  my $q_vec_file = "";
  if(exists $cmnd_hashref->{"q-vec-file"}) {
    $q_vec_file = $cmnd_hashref->{"q-vec-file"};
  }
  my $x_guess_file = "";
  if(exists $cmnd_hashref->{"x-guess-file"}) {
    $x_guess_file = $cmnd_hashref->{"x-guess-file"};
  }
  my $p_guess_file = "";
  if(exists $cmnd_hashref->{"p-guess-file"}) {
    $p_guess_file = $cmnd_hashref->{"p-guess-file"};
  }
  my $scale_x_guess = 1.0;
  if(exists $cmnd_hashref->{"scale-x-guess"}) {
    $scale_x_guess = $cmnd_hashref->{"scale-x-guess"};
  }
  my $scale_p_guess = 1.0;
  if(exists $cmnd_hashref->{"scale-p-guess"}) {
    $scale_p_guess = $cmnd_hashref->{"scale-p-guess"};
  }
  my $x_solu_file = "";
  if(exists $cmnd_hashref->{"x-solu-file"}) {
    $x_solu_file = $cmnd_hashref->{"x-solu-file"};
  }
  my $p_solu_file = "";
  if(exists $cmnd_hashref->{"p-solu-file"}) {
    $p_solu_file = $cmnd_hashref->{"p-solu-file"};
  }
  my $black_box = (exists $cmnd_hashref->{"black-box"});
  my $fd_step_len = -1.0;
  if(exists $cmnd_hashref->{"fd-step-len"}) {
    $fd_step_len = $cmnd_hashref->{"fd-step-len"};
  }
  #
  return
    ($fd_step_len > 0.0 && $use_finite_diff_obj && $use_finite_diff_con ? " --no-support-derivatives" : "" )
    .($q_vec_file ne "" ? " --q-vec-file=$q_vec_file" : "" )
    ." --extra-moocho-thyra-params=\""
       ."<ParameterList>"
         .($do_sim ? "<Parameter name=\\\"Solve Mode\\\" type=\\\"string\\\" value=\\\"Forward Solve\\\"/>" : "" )
         .($use_direct ? "<Parameter name=\\\"NLP Type\\\" type=\\\"string\\\" value=\\\"Direct\\\"/>" : "" )
         .(!$do_sim && $use_default_inv_obj
           ? "<Parameter name=\\\"Use Built-in Inverse Objective Function\\\" type=\\\"bool\\\" value=\\\"1\\\"/>"
             ."<ParameterList name=\\\"Inverse Objective Function Settings\\\">"
               ."<Parameter name=\\\"Observation Multiplier\\\" type=\\\"double\\\" value=\\\"1.0\\\"/>"
               ."<ParameterList name=\\\"Observation Target Vector\\\">"
                 ."<Parameter name=\\\"File Name Base\\\" type=\\\"string\\\" value=\\\"$q_vec_file\\\"/>"
               ."</ParameterList>"
               ."<Parameter name=\\\"Parameter Multiplier\\\" type=\\\"double\\\" value=\\\"0.0\\\"/>"
               ."<ParameterList name=\\\"Parameter Base Vector\\\">"
                 ."<Parameter name=\\\"File Name Base\\\" type=\\\"string\\\" value=\\\"$p_guess_file\\\"/>"
               ."</ParameterList>"
             ."</ParameterList>"
           : ""
           )
         .($x_guess_file ne ""
           ? "<ParameterList name=\\\"State Guess\\\">"
               ."<Parameter name=\\\"File Name Base\\\" type=\\\"string\\\" value=\\\"$x_guess_file\\\"/>"
               .($scale_x_guess != 1.0 ? "<Parameter name=\\\"Scale By\\\" type=\\\"double\\\" value=\\\"$scale_x_guess\\\"/>" : "" )
             ."</ParameterList>"
           : ""
           )
         .($p_guess_file ne ""
           ? "<ParameterList name=\\\"Parameter Guess\\\">"
               ."<Parameter name=\\\"File Name Base\\\" type=\\\"string\\\" value=\\\"$p_guess_file\\\"/>"
               .($scale_p_guess != 1.0 ? "<Parameter name=\\\"Scale By\\\" type=\\\"double\\\" value=\\\"$scale_p_guess\\\"/>" : "" )
             ."</ParameterList>"
           : ""
           )
         .($x_solu_file ne "" ? "<Parameter name=\\\"State Solution File Base Name\\\" type=\\\"string\\\" value=\\\"$x_solu_file\\\"/>" : "" )
         .($p_solu_file ne "" ? "<Parameter name=\\\"Parameters Solution File Base Name\\\" type=\\\"string\\\" value=\\\"$p_solu_file\\\"/>" : "" )
         .($black_box ? "<Parameter name=\\\"Nonlinearly Eliminate States\\\" type=\\\"bool\\\" value=\\\"1\\\"/>" : "" )
         .($black_box ? "<Parameter name=\\\"Forward Newton Tolerance\\\" type=\\\"double\\\" value=\\\"$fwd_newton_tol\\\"/>" : "" )
         .($fd_step_len > 0.0
           ?
           ($use_finite_diff_obj
            ? (
               "<Parameter name=\\\"Use Finite Differences For Objective\\\" type=\\\"bool\\\" value=\\\"1\\\"/>"
               ."<ParameterList name=\\\"Objective Finite Difference Settings\\\">"
               ."<Parameter name=\\\"FD Step Length\\\" type=\\\"double\\\" value=\\\"$fd_step_len\\\"/>"
               ."</ParameterList>"
               )
            : ""
            )
           .($use_finite_diff_con
             ? (
                "<Parameter name=\\\"Use Finite Differences For Constraints\\\" type=\\\"bool\\\" value=\\\"1\\\"/>"
                ."<ParameterList name=\\\"Constraints Finite Difference Settings\\\">"
                ."<Parameter name=\\\"FD Step Length\\\" type=\\\"double\\\" value=\\\"$fd_step_len\\\"/>"
                ."</ParameterList>"
                )
             : ""
             )
           : ""
           )
       ."</ParameterList>\"";
}

sub run_cmnd {
  my $cmnd = shift;
  $cmnd .= " 2>&1 | tee run-test.out";
  print "\nRunning command: $cmnd\n";
  my $rtn = system($cmnd);
  my $passed_msg = `grep 'End Result: TEST PASSED' run-test.out`;
  chomp($passed_msg);
  print "\npassed_msg = \"$passed_msg\"\n";
  return ( $passed_msg eq "End Result: TEST PASSED" ? 0 : 1);
}

sub run_test {
  my $cmnd = shift;
  my $curr_dir = cwd();
  print
    "\n***"
    ,"\n*** $curr_dir"
    ,"\n***\n";
  return run_cmnd($cmnd);
}

sub run_case {

  my $cmnd = shift;
  my $do_fwd_prob = shift;
  my $do_inv_prob = shift;
  my $fwd_dir = shift;

  $do_fwd_prob = 1 if(!defined($do_fwd_prob));
  $do_inv_prob = 1 if(!defined($do_inv_prob));
  $fwd_dir = cwd()."/fwd" if (!defined($fwd_dir));

  my $rtn = 0;

  if($solve_fwd_prob) {

    if($do_fwd_prob) {

      mkchdir("fwd-init");

      my $fwdcmnd =
        $cmnd
        ." --moocho-options-file=$fwd_init_options_file "
        .moocho_thyra_solver_cmndline(
            {
              "do-sim"=>""
              ,"use-first-order"=>""
              ,"x-solu-file"=>"x.out"
              ,"p-solu-file"=>"p.out"
            }
          );
      $rtn = 1 if (0!=run_test($fwdcmnd));

      chdir "..";

      mkchdir("fwd");

      $fwdcmnd =
        $cmnd
        ." --moocho-options-file=$fwd_options_file "
        .moocho_thyra_solver_cmndline(
            {
              "do-sim"=>""
              ,"use-first-order"=>""
              ,"x-guess-file"=>"../fwd-init/x.out"
              ,"p-guess-file"=>"../fwd-init/p.out"
              ,"scale-p-guess"=>$p_solu_scale
              ,"x-solu-file"=>"x.out"
              ,"p-solu-file"=>"p.out"
            }
          );

      $rtn = 1 if (0!=run_test($fwdcmnd));

      chdir "..";

    }

    if($do_inv_prob) {

      my $dir_name = "inv";
      mkdir $dir_name, 0777;
      chdir $dir_name;

      my $build_subdirs =
        (
         ( $do_analytic ? 1 : 0 )
         +
         ( $do_black_box ? 1 : 0 )
         +
         ( $do_finite_diff ? 1 : 0 )
        ) > 1;

      if($do_black_box) {
        mkchdir("black-box") if($build_subdirs);
        my $invcmnd =
          $cmnd 
            ." --moocho-options-file=$bb_options_file "
            .moocho_thyra_solver_cmndline(
               {
                 "q-vec-file"=>"$fwd_dir/../fwd-init/x.out"
                 ,"black-box"=>""
                 ,"fd-step-len"=>$bb_fd_step_len
                 ,"x-guess-file"=>"$fwd_dir/x.out"
                 ,"p-guess-file"=>"$fwd_dir/p.out"
               }
             );
        $rtn = 1 if (0!=run_test($invcmnd));
        chdir("..") if($build_subdirs);
      }

      if($do_finite_diff) {
        mkchdir("finite-diff") if($build_subdirs);
        my $invcmnd =
          $cmnd 
            ." --moocho-options-file=$fd_options_file "
            .moocho_thyra_solver_cmndline(
               {
                 "q-vec-file"=>"$fwd_dir/../fwd-init/x.out"
                 ,"fd-step-len"=>$fd_step_len
                 ,"x-guess-file"=>"$fwd_dir/x.out"
                 ,"p-guess-file"=>"$fwd_dir/p.out"
               }
             );
        $rtn = 1 if (0!=run_test($invcmnd));
        chdir("..") if($build_subdirs);
      }

      if($do_analytic) {
        mkchdir("analytic") if($build_subdirs);
        my $invcmnd =
          $cmnd 

            ." --moocho-options-file=$inv_options_file "
            .moocho_thyra_solver_cmndline(
               {
                 "q-vec-file"=>"$fwd_dir/../fwd-init/x.out"
                 ,"x-guess-file"=>"$fwd_dir/x.out"
                 ,"p-guess-file"=>"$fwd_dir/p.out"
               }
             );
        $rtn = 1 if (0!=run_test($invcmnd));
        chdir("..") if($build_subdirs);
      }

      chdir "..";

    }

  }
  else {

    $rtn = 1 if (0!=run_test($cmnd));

  }

  return $rtn;

}
