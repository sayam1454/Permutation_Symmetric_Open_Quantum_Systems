
# Introduction
This work aims to simulate open quantum systems for an arbitary d-level quantum system.The task to simulate the dynamics of a quantum system when interacting with its environment is numerically inefficient. The scaling with system size grows exponentialy (eg: For qubits: O(4^N)). 
<br>
To overcome this challenge, one can invoke a particular symmetry called the permutation symmetry. This assumption brings down the scaling from exponential to polynomial (eg: For qubits: O(N^3) and in some special case O(N^2)). 
<br>
Not only it is mathematically helpful, this kind of formalism can be used to describe various experimental setups that are of major interests.(eg: all to all coupled systems)

<section>
  <h1>Permutation-symmetric dissipative processes</h1>
  <p>
    Dissipative processes in open quantum systems can be divided into
    <strong>local</strong> and <strong>collective</strong> categories. 
    Local processes act on each particle or site independently, while 
    collective processes act on the entire ensemble through shared
    interactions with a common environment.
  </p>

  <h2>Local processes</h2>
  <p>
    Local processes describe dynamics where each atom, spin, or qutrit
    interacts with its own environment. They generally break permutation 
    symmetry because each particle can evolve differently, even if the 
    rates are the same.
  </p>
  <ul>
    <li>
      <strong>Spontaneous decay:</strong> Each site loses its excitation
      independently by emitting energy into its own bath. This reduces the
      excited population but does not create collective effects such as
      superradiance.
    </li>
    <li>
      <strong>Pumping:</strong> Each site can be excited individually by
      absorbing energy from a local bath or noise source. The excitation
      happens independently across the ensemble.
    </li>
    <li>
      <strong>Decoherence (dephasing):</strong> Random phase noise acts
      independently on each site, destroying local quantum coherences and
      reducing the possibility of collective interference effects.
    </li>
  </ul>

  <h2>Collective processes</h2>
  <p>
    Collective processes occur when all sites couple to the same bath or 
    field. These processes preserve permutation symmetry and act within the 
    symmetric (Dicke-like) subspace. As a result, they can generate 
    genuinely collective phenomena such as superradiance or subradiance.
  </p>
  <ul>
    <li>
      <strong>Collective decay:</strong> Excitations are released into a 
      shared environment. When multiple excitations are present, the decay
      rate can become enhanced (superradiance) or suppressed (subradiance).
    </li>
    <li>
      <strong>Collective pumping:</strong> The entire ensemble is excited 
      in a correlated way through a common bath or shared drive. This leads 
      to collective population of excited states rather than independent 
      excitations.
    </li>
    <li>
      <strong>Collective dephasing:</strong> All sites experience the same 
      phase noise, such as from a fluctuating global field. This preserves 
      symmetry because the relative phases between sites remain unaffected, 
      although coherences between different energy sectors are suppressed.
    </li>
  </ul>

  <h2>Rates in permutation-symmetric systems</h2>
  <p>
    In permutation-symmetric systems, all particles are identical and
    indistinguishable. For this reason:
  </p>
  <ul>
    <li>
      The <strong>local rates</strong> (decay, pump, dephasing) are the same
      for each site, since every particle experiences the same environment.
    </li>
    <li>
      The <strong>collective rates</strong> (decay, pump, dephasing) apply
      equally to the entire ensemble, as all particles couple to the same
      bath or field.
    </li>
  </ul>
  <p>
    This equality of rates ensures that the dynamics remain permutation 
    symmetric: the system can be described without distinguishing between 
    individual sites, which greatly simplifies simulations in the symmetric 
    subspace.
  </p>

# Benchmarking
The intractable nature of this problem and the specific symmetry makes this a relatively less explored area. However ther is an existing framework in python incorporaated with QuTip called PIQS (Permutation Symmetry Quantum Solver). There is however no efficient framework for qudit system in python.
<br>
The method we employ is adopted from Gegg, M. (2017). Identical emitters, collective effects and dissipation in quantum optics: Novel numerical approaches for quantum master equations (Order No. 27610221). Available from ProQuest Dissertations & Theses Global. (2424459914). doi:https://doi.org/10.14279/depositonce-6526 Retrieved from https://www.proquest.com/dissertations-theses/identical-emitters-collective-effects-dissipation/docview/2424459914/se-2. Necessary fine tuning are made whenever felt so.
<br>
The first benchmarking is done by comparing this code for qubits with the performance of PIQS. 
<br>
For qudits however we dont have a dedicated python package. So for that we work in completely uncoupled basis in QuTip. Though it is inefficient, this allows us to compare the code for a very small sample size and gain maximum confidence in this method.
<br>
<br>
<strong>Note</strong>: There exist a C++ library mady by the author of the above mentioned paper called PsiQuaSP. However thepopularity of python for variety of nature, calls for a python based simulation code for such permutation symmetric d-level states.