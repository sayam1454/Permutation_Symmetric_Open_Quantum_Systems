# Introduction
This work aims to simulate open quantum systems for an arbitary d-level quantum system.The task to simulate the dynamics of a quantum system when interacting with its environment is a numerically inefficient task. The scaling with system size grows as O(4^N). 
<br>
To overcome this challenge, one can invoke a particular symmetry called the permutation symmetry. This assumption brings down the scaling from exponential to polynomial O(N^3) and in some special case O(N^2). 

The code models for both the local and the collective dissipative terms. A  small explaination of all these terms and how they behave in a permutation symmetric system are as follows:
<br>
<ol>
<li><b>LOCAL PROCESSSES </b></li>

Local processes are
<ol>
<li><b>Spontaneous Decay </b>: </li>
<li><b>Pumping</b>: </li>
<li><b>Decoherence</b>: </li>
</ol>

<br>

<li><b>COLLECTIVE PROCESSSES </b></li>

Collective processes are
<ol>
<li><b>Collective Decay </b>: </li>
<li><b>Collective Pumping </b>: </li>
<li><b>Collective Decoherence </b>: </li>
</ol>
</ol>

# Benchmarking
The existing PIQS package QuTip is the currently widely used package to simulate permutationally symmetric open quantum dynamics of 