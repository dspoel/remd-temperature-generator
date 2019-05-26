<?php
include "functions.inc";
showTitle();

// Numeric fields
define("Pdes",$_POST['Pdes']);
define("Tlow",$_POST['Tlow']);
define("Thigh",$_POST['Thigh']);
define("Nw",$_POST['Nw']);
define("Np",$_POST['Np']);
define("Tol",$_POST['Tol']);

// Multiple choice fields
define("PC",$_POST['PC']);
define("WC",$_POST['WC']);
define("Hff",$_POST['Hff']);
define("Vs",$_POST['Vs']);
define("Alg",$_POST['Alg']);

// Constants. May depend on input in principle (but not yet).
define("A0","-59.2194");
define("A1","0.07594");
define("B0","-22.8396");
define("B1","0.01347");
define("D0","1.1677");
define("D1","0.002976");
define("maxiter",100);

$kB       = 0.008314;
$Tlow     = Tlow;
$Thigh    = Thigh;
$Nw       = Nw;
$Np       = Np;
$Npp      = 0;
$Hff      = Hff;
$Vs       = Vs;
$Pdes     = Pdes;
$PC       = PC;
$NC       = 0;
$VC       = 0;
$Tol      = Tol;
$debug    = 0;

// Check input variables
$error = 0;
if ( !is_numeric($Pdes) || !is_numeric($Tlow) || !is_numeric($Thigh) ||
     !is_numeric($Np) || !is_numeric($Nw) || !is_numeric($Tol) ) {
  print "<p><b>ERROR</b>: Some of your inputs are not numbers</p>\n";
  $error++;
}
if ( ($Pdes > 1) || ($Pdes < 0) ) {
  print "<p><b>ERROR</b>: You have to give a probability P<sub>des</sub>: between 0 and 1!</p>\n";
  $error++;
}
if ( Thigh <= Tlow ) {
  print "<p><b>ERROR</b>: The lower limit of the temperature range has to be below the upper limit - Check your input!</p>\n";
  $error++;
}
if ( (Tlow <= 0) || (Thigh <= 0) ) {
  print "<p><b>ERROR</b>: You must have temperatures that are &gt; 0!</p>\n";
  $error++;
}
if ( $Np==0 ) {
  print "<p><b>ERROR</b>: You can not have zero atoms in protein!</p>\n";
  $error++;
}
if ( Alg != 0) {
  print "<p><b>ERROR</b>: Can not do constant volume yet!</p>\n";
  $error++;
}

function calc_mu($Nw,$Np,$Temp,$FEner) {
  return ((A0+A1*$Temp)*$Nw + 
	  (B0+B1*$Temp)*$Np -
	  $Temp*$FEner);
}

function erf($x) {
    # constants
    $a1 =  0.254829592;
    $a2 = -0.284496736;
    $a3 =  1.421413741;
    $a4 = -1.453152027;
    $a5 =  1.061405429;
    $p  =  0.3275911;

    # Save the sign of x
    $sign = 1;
    if ($x < 0) {
        $sign = -1;
    }
    $x = abs($x);

    # Abramowitz, M. and Stegun, I. A. 1972. 
    # Handbook of Mathematical Functions with Formulas, Graphs, and Mathematical Tables. Dover.
    # Formula 7.1.26
    $t = 1.0/(1.0 + $p*$x);
    $y = 1.0 - ((((($a5*$t + $a4)*$t) + $a3)*$t + $a2)*$t + $a1)*$t*exp(-$x*$x);

    return $sign*$y;
}

function myeval($m12, $s12, $CC, $u) {
  $argument = -$CC*$u - ($u-$m12)**2/(2*$s12**2);
  return exp($argument);
}

function myintegral($m12, $s12, $CC) {
  $int  = 0.0;
  $umax = $m12+5*$s12;
  $du   = 0.5;
  if ($debug > 1) {
    print("umax = $umax m12 = $m12 s12 = $s12 CC = $CC<br>");
  }
  for($u = 0; $u<$umax; $u+=$du) {
    $di = myeval($m12, $s12, $CC, $u+$du/2);
    $int += $di;
  }
  $pi = 3.14159265358979;
  return $du*$int/($s12*sqrt(2*$pi));
}

if ($error == 0) {
  $Npp = 0;
  $Nprot = 0;
  if ($Hff == 0) {
    $Nh = round($Np*0.5134);
    if (Vs == 1) {
      $VC = round(1.91*$Nh);
    }
    $Nprot = $Np;
  }
  else {
    $Npp = round($Np/0.65957);
    $Nh = round($Np*0.22);
    if (Vs == 1) {
      $VC = round($Np+1.91*$Nh);
    }
    $Nprot = $Npp;
  }
  if ($PC == 1) {
    $NC = $Nh;
  }
  else if ($PC == 2) {
    $NC = $Np;
  }
  
  $Ndf      = (9-WC)*$Nw + 3*$Np-$NC-$VC;
  $FlexEner = 0.5*$kB*($NC+$VC+WC*$Nw);

  /////////////////////
  print <<<_LAYOUT1_
    <html>
    <p>
    <h2>Summary of input and derived variables.</h2>
    <table border=2 cellpadding=3 bordercolor=yellowgreen>
    <tr><th>Variable</th><th>Value</th></tr>
    <tr><td>P<sub>des</sub></td><td>$Pdes</td></tr>
    <tr><td>Temperature range</td><td>$Tlow - $Thigh</td></tr>
    <tr><td>Number of water molecules</td><td>$Nw</td></tr>
    <tr><td>Number of protein atoms</td><td>$Np</td></tr>
_LAYOUT1_;

  if ( $Npp > 0) {
    print <<< _LAYOUT2_
      <tr><td>Including all H</td><td>~ $Npp</td></tr>
_LAYOUT2_;

  }
  print <<< _LAYOUT3_
    <tr><td>Number of hydrogens in protein</td><td>~ $Nh</td></tr>
    <tr><td>Number of constraints</td><td>~ $NC</td></tr>
    <tr><td>Number of vsites</td><td>~ $VC</td></tr>
    <tr><td>Number of degrees of freedom</td><td>~ $Ndf</td></tr>
_LAYOUT3_;
  printf ("<tr><td>Energy loss due to constraints</td><td>%.2f (kJ/mol K)</td></tr> </table> </p>\n",$FlexEner);

/////////////////////

  $index = 1;
  $T[$index] = Tlow;
  
  while ( ($T[$index] < Thigh) ) {
    $piter   = 0;
    $forward = 1;
    $iter    = 0;
    $T1      = $T[$index];
    $T2      = $T1+5;
    if ( $T2 >= Thigh ) { 
       $T2 = Thigh;
    }
    $low     = $T1;
    $high    = Thigh;
    if ($debug == 2) {
      printf("<p>Index %d, T1 = %f, T2 = %f</p>\n", $index, $T1, $T2);
    }
    while ( (abs($Pdes-$piter) > $Tol) && ($iter < maxiter) ) {
      $iter++;
      $mu12 = ($T2-$T1) * ((A1*$Nw)+(B1*$Nprot)-$FlexEner);
      $MM[$index] = $mu12;
      
      $CC    = (1/$kB) * ( (1/$T1)-(1/$T2) );
      $Delta = $CC*$mu12;
      
      $var = $Ndf*(D1*D1*( $T1*$T1 + $T2*$T2 ) +
		   2*D1*D0*($T1+$T2) +
		   2*D0*D0);
      
      $sig12 = sqrt($var);
      $SS[$index] = $sig12;
      
      if ($sig12 == 0) {
	print "Sigma = 0\n";
	exit(1);
      }
      // I1
      $erfarg1 = -$mu12/($sig12*sqrt(2));
      $I1      = 0.5*(1 + erf($erfarg1));

      // I2
      // Old analytical code according to the paper, however
      // this suffers from numerical issues in extreme cases.
      // $exparg  = $CC*(-$mu12 + $CC*$var/2);
      // $erfarg2 = ($mu12 - $CC*$var)/($sig12*sqrt(2));
      // $I2      = 0.5*exp($exparg)*(1.0 + erf($erfarg2));
      // Use numerical integration instead.
      $I2      = myintegral($mu12, $sig12, $CC);
      $piter   = ($I1 + $I2);
      
      if ($debug == 2) {
         printf("<p>\n");
      	printf("DT = %.3f CC = %.4f<br>", $T2-$T1, $CC);
	printf("mu12 = %.1f, sig12 = %.1f, Delta = %.1f<br>",
               $mu12, $sig12, $Delta);
	printf("erfarg1 = %.4f, Integral1 = %.4f<br>", $erfarg1, $I1);
	//printf("exparg = %.2e, erfarg2 = %.1f<br>", $exparg, $erfarg2);
        //printf("myintegral = %.4f, Integral2 = %.4f<br>", $myint, $I2);
	printf("piter = %.3f<br>", $piter);
        printf("</p>\n");
      }
      
      if ( $piter > $Pdes ) {
	if ( $forward==1 ) {
	  $T2 = $T2 + 1.0; 
	}
	elseif ( $forward==0 ) {
	  $low = $T2;
	  $T2 = $low + (($high-$low)/2);
	}
        if ( $T2 >= Thigh ) { 
           $T2 = Thigh;
        }      
      }
      elseif ( $piter < $Pdes ) {
	if ( $forward==1 ) {
	  $forward = 0;
	  $low = $T2 - 1.0;
	}
	$high = $T2;
	$T2 = $low + (($high-$low)/2);
      }
    }
    
    $P[$index]      = $piter;
    $Siigma[$index] = sqrt($Ndf)* (D0 + D1*$T1);
    $Muu[$index]    = calc_mu($Nw,$Nprot,$T1,$FlexEner);
    
    $index++;
    $T[$index] = $T2;
  }
  $Siigma[$index] = sqrt($Ndf)* (D0 + D1*$T[$index]);
  $Muu[$index] = calc_mu($Nw,$Nprot,$T[$index],$FlexEner);


/////////////////////
  print <<<_YOGHURT_
    <hr noshade>
    <h2>Temperatures and Energies.</h2>
    <P>
    <table border=2 cellpadding=3>
    <tr>
    <th></th><th>Temperature (K)</th><th>&mu; (kJ/mol)</th><th>&sigma; (kJ/mol)</th><th>&mu;<sub>12</sub> (kJ/mol)</th><th>&sigma;<sub>12</sub> (kJ/mol)</th><th>P<sub>12</sub></th>
										      </tr>
_YOGHURT_;
    for($k=1; ($k<=$index); $k++) {
      printf("<tr><td align=center>%d</td><td>%.2f</td>",$k,$T[$k]);
      if ($k == 1)
	printf("<td>%.0f</td><td>%.2f</td><td></td><td></td><td></td>",
	       $Muu[$k],$Siigma[$k]);
      else {
	printf("<td>%.0f</td><td>%.2f</td><td>%.1f</td><td>%.2f</td><td>%.4f</td>",
	       $Muu[$k],$Siigma[$k],$MM[$k-1],
	       $SS[$k-1],$P[$k-1]);
      }
      printf("</tr>\n");
    }  
    
    print <<<_LAYOUT2_
      </table>
      </P>
      <hr noshade>
      </body>
      </html>
_LAYOUT2_;
/////////////////////
    printf("<p>For your scripting pleasures we also give the temperatures below as one long comma-separated string. Enjoy.<br></p>\n<p>\n");
    for($k=1; ($k<$index); $k++) {
      printf("  %.2f,",$T[$k]);
    }
    printf("  %.2f\n</p>\n",$T[$k]);
}

showFooter();

?>
