<?php

function showTitle()
{
print <<<_PONY_
<html>
<head><title>Temperature generator for REMD-simulations</title></head>
<font color=darkgreen>
_PONY_;
}

function showIntro()
{
print <<<_PAARD_
<FONT COLOR = "darkgreen"> 
<H1>Temperature generator for REMD-simulations</H1>

<P>
This site contains a webserver for generating temperatures for
REMD-calculations. You submit the number of protein atoms and water
molecules in your system, and an upper and lower limit for the
temperature range, information about constraints and/or virtual sites
and a desired exchange probability P<sub>des</sub>, and the webserver
will predict a temperature series with correspondig energy differences
and standard deviations which matches the desired probability
P<sub>des</sub>.  You can then use these temperatures in REMD
simulations.
</P>

<p> A word of caution is in place here. The derivation of the
parameters for the prediction was done with the OPLS/AA force field
and the 
<a href="http://www.gromacs.org">GROMACS</a>
 software. When
using other force fields, software and/or other algorithms (cut-off
treatment, pressure and temperature scaling etc.) results may deviate,
although the tests performed in our paper, including using the GROMOS96
force field, show that these deviations are minor. Nevertheless you
are encouraged to check your exchange probabilities, and compare them
to the desired probabilities.  </p>

<hr noshade>
_PAARD_;
}

function showForm()
{
print <<<_KOE_
<form action="tgenerator.php" method="post">
<table border="2" bordercolor="darkgreen">

<tr>
<td>Exchange probability:</td>
<td><input type="text" name="Pdes"></td>
<td>Tolerance:</td>
<td><input type="text" name="Tol" value="1e-4"></td>
</tr>

<tr>
<td>Lower temperature limit:</td>
<td> <input type="text" name="Tlow"></td>
<td>Upper temperature limit:</td>
<td> <input TYPE="TEXT" NAME="Thigh"></td>
</tr>

<tr>
<td>Number of water molecules:</td>
<td><input TYPE="TEXT" NAME="Nw"></td>
<td>Constraints in water:</td>
<td> <SELECT NAME="WC">
<OPTION VALUE="0">Fully Flexible
<OPTION VALUE="2">Flexible angle
<OPTION VALUE="3">Rigid
</SELECT></td>
</tr>

<tr>
<td>Number of protein atoms:</td>
<td> <input TYPE="TEXT" NAME="Np"></td>
<td>Constraints in the protein:</td>
<td> <SELECT NAME="PC">
<OPTION VALUE="0">Fully Flexible
<OPTION VALUE="1">Bonds to hydrogens only
<OPTION VALUE="2">All bonds
</SELECT></td>
</tr>

<tr>
<td>Hydrogens in protein:</td>
<td>  <SELECT NAME="Hff">
<OPTION VALUE="0">All H
<OPTION VALUE="1">Polar H
</SELECT></td>
<td>Virtual sites in protein:</td>
<td>  <SELECT NAME="Vs">
<OPTION VALUE="0">None
<OPTION VALUE="1">Virtual Hydrogen
</SELECT></td>
</tr>

<tr>
<td>Simulation type:</td>
<td>  <SELECT NAME="Alg">
<OPTION VALUE="0">NPT
<OPTION VALUE="1">NVT
</SELECT></td>
<td></td><td></td>
</tr>
</table>

<p>
<INPUT TYPE="SUBMIT" VALUE="Submit">
<INPUT TYPE="RESET" VALUE="Clear">
</p>

</form>
_KOE_;
}

function showCounter()
{
  print "<IMG SRC=\"http://folding.bmc.uu.se/cgi-bin/nph-count?width=6&link=index.html\">";
  print "<HR NOSHADE>\n";
}
function showFooter()
{
print <<<_SCHAAP_
<HR NOSHADE>

<p>
If you use the results from this webserver in simulations which are published in scientific journals, please
cite:<br> Alexandra Patriksson and David van der Spoel, <i>A
temperature predictor for parallel tempering simulations</i>
Phys. Chem. Chem. Phys.,  10 pp. 2073-2077 (2008) <a href="http://dx.doi.org/10.1039/b716554d">http://dx.doi.org/10.1039/b716554d</a>.
</p>
<p>

We also recommend the following literature about theory behind replica
exchange simulations [1,2] and applications of REMD [3,4]. A recent
review about sampling is in ref. [5].

<ol>

<li>K. Hukushima and K. Nemoto: <i>Exchange Monte Carlo Method and
Application to Spin Glass Simulations</i> J. Phys. Soc. Jpn. <b>65</b>
pp. 1604-1608 (1996)</li>

<li>T. Okabe and M. Kawata and Y. Okamoto and M. Mikami:
<i>Replica-exchange {M}onte {C}arlo method for the isobaric-isothe
rmal ensemble</i> Chem. Phys. Lett. <b>335</b> pp. 435-439 (2001)</li>

<li>Marvin Seibert, Alexandra Patriksson, Berk Hess and David van der
Spoel: <i>Reproducible polypeptide folding and structure prediction
using molecular dynamics simulations</i> J. Mol. Biol. <b>354</b>
pp. 173-183 (2005)</li>

<li>David van der Spoel and M. Marvin Seibert: <i>Protein Folding
Kinetics and Thermodynamics from Atomistic Simulations</i>
Phys. Rev. Lett. <b>96</b> pp. 238102 (2006)</li>

<li>H. X. Lei and Y. Duan: <i>Improved sampling methods for molecular
simulation</i> Curr. Opin. Struct. Biol. <b>17</b> pp. 187-191
(2007)</li>

</ol>
</p>
<p>
In case of questions please mail to <b>t-remd</b> at <b>xray.bmc.uu.se</b>.
</p>
</font>
<hr noshade>
</html>
_SCHAAP_;
}

?>
