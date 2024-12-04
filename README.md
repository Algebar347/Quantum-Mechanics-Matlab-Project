# Quantum-Mechanics-Matlab-Project
#### Computational project for PC4230 Quantum MechanicsⅢ. Including topics: DVR, FFT, perturbation theory, transition probability, one-photon process, et al.

<br/>

>Consider a particle moving in a harmonic oscillator potential, with the associated dimensionless Hamiltonian given by:
>
>$` H^0 = \frac{p^2}{2}+\frac{1}{2} x^2 `$  with the effective  $`\hbar = 1`$. 
>
>In terms of the dimensionless units adopted here, the energy level spacing is 1.0 and initially the system is in its ground state. We now switch on a timeperiodic perturbation $`V(t)`$ (t is in dimensionless units here), whose dimensionless form is given by:
>
>$`V(t) = A \sin{x} \cos{\omega t}`$, with A being the driving amplitude and $`\omega`$ the driving frequency. 
>
>To examine the ensuing dynamics you are now asked to adapt the code you learned before to implement the split-operator method for time-dependent Hamiltonians. To that end you need to divide the time into many small intervals, and within each small time window, you can still approximate the Hamiltonian as a constant and then use the split-operator method to obtain the evolution operator for each small time interval.

1. Please adapt the circulated codes on Canvas to **simulate the time evolution of this system**. Before you do any physical investigations, make sure that your results have convergence (that is, physics observed is no longer dependent on time step size, number  of DVR points etc). **Make interesting observations by tuning A and ω**.

2. You may use analytical forms of the eigenstates of a harmonic oscillator available from internet (note that Hermite polynomials can be directly called in Matlab) when you need to define the initial state or when you analyze the quantum amplitude projected onto the eigenstates of the unperturbed Harmonic oscillator system. You may also numerically **generate the ground state and all other excited states from the DVR approach we studied**. For the second approach, some caution is needed as the convention (and boundary conditions) in our lectures for **defining the grid points** used in DVR for spectrum calculations and that for FFT calculations are slightly different and so you need to do something about it.

3. Please **analyze the transition probabilities** from the ground state to the first excited state when the driving frequency is almost on-resonance with the natural frequency of the harmonic oscillator. In particular, how does the transition probability depend on time and on the driving frequency? You are encouraged to discuss any feature that is interesting to you. Remember, you are actually doing an experiment and so feel free to explore.

4. Compare your results from step 3 with our **first-order time-dependent perturbation theory** with or without the so-called rotating-wave approximation. Analyze how firstorder time-dependent perturbation theory breaks down as A increases.

5. As we increase the driving field amplitude A, is there any **two-photon transition probability** from the ground state to the first excited state if ω ≈ 0.5? That is, are the transitions from the ground state to the first excited state due to off-resonance onephoton processes or due to any two-photon processes? What is your reasoning?

6. Consider now ω ≈ 2.0, which is almost equal to the transition frequency between the  ground state and the second excited state. What would the first-order perturbation theory predict regarding the transition probabilities **from the ground state to the second excited state**? What do you see from your numerical simulations? What about the transition probabilities to the fourth excited state when ω ≈ 2.0 – How to explain your simulation results in terms of some theory we studied already?

7. (More adventure). You are encouraged to go beyond what is assigned above (e.g., introducing nonlinearities to the oscillator, considering multiple driving fields, and considering a very small driving frequency etc). Remember, what you have now is basically a small computational lab and how much you can “discover” is mainly limited  by how much you are willing to try.
