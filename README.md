This repository contains three scripts.

Phase_demodulation_skript_complete: This is meant to post-process the experimentally derived interferograms, recorded during the measurement campaign. It is able to demodulate the phase and subsequently unwrap it. The demodulation algorithm of choice is the windowed Fourier transforn (WFT). The phase unwrapping is performed by a weighted phase unwrapping algorithm. It is therefore important, that the script 'weighted_phase_unwrap' is added to the MatLab data path. Afterwards the phase differences are calculated from the gradients via integration. The resulting phase field is then transformed into a temperature field by scaling it to the experimental boundary conditions. This script yields:

masked versions of the input interferograms for qualitative analysis
demodulated and unwrapped phase of said interferograms
the resulting tempersture profile, which has been calculated from the demodulated phase
Two exemplary interferograms have been attached ('2K_Hor_X2', '2K_Hor_Y2').

Simple_phase_demodulation: A generalised and simplified version of the previous script, which can be used to produce a simple, wavy, synthetic interferogram. This interferogram is then demodulated, using the WFT algorithm. Afterwards a phase unwrap is performed. Again the script 'weighted_phase_unwrap' needs to be added to the MatLab data path.

Artificial_interferogram_generator: This script is used to calculate circular synthetic temperature profile from an input 3D temperature field by integrating along the z-axis. This temperature profile is then modulated into a set of synthetic interferograms for the x- and y-gradient respectively.

An exemplary 3D temperature field can be found here: https://www.b-tu.de/owncloud/s/X9nzJpDHgdjqkMx
