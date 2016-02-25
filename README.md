# SpectraJu

Copyright (c) 2016 Charles Le Losq

email: charles.lelosq@anu.edu.au

Licence MIT: see LICENCE.md

SpectraJu is a package aimed at helping spectroscopic (Raman, Infrared, Nuclear Magnetic Resonance, XAS...) data treatment under Julia.

It's aim is to provide the simplest way to perform actions like baseline fitting and removal or peak fitting for instance, while respecting the freedom offered by data treatment through coding.

It is particularly focused on large datasets because of the high speed of Julia's, e.g. for performing peak fitting along Infrared diffusion profiles.

For peak fitting, the JuMP interface offers a very flexible yet clear way to build models, that can be solve with top-notch solvers such as Ipopt.

Examples are included as notebooks.

Package is under construction, any help welcome!

Last modified 26/02/2016




