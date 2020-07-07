<script src="https://polyfill.io/v3/polyfill.min.js?features=es6"></script>
<script id="MathJax-script" async src="https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-mml-chtml.js"></script>
# Ion optics

### Optimization of Twiss parameters to influence the emissivity of an ion beam

<p>Programs developed during an internship in ion optics to discover the notion of emissivity, Twiss parameters and "sigma" matrices for ion beams passing through different optical systems (drift, thin lens, Einzel lens, magnetic dipole, electrostatic quadrupole).</p>

<p>The parameters are given in **parameters.py**, which allows to change the values of the optical elements. The Twiss parameter "gamma" is calculated in **gamma.py**.</p>

$$
\gamma={\frac {\epsilon^2+\alpha^2 \epsilon^2}{\beta \epsilon^2}}
\tag{1}
\label{1}
$$

<p>The "trash cleanup" part allows to remove (or not) files that have been in the trash for more than 1 week. But you can change this value with the corresponding macro. Explanations here: <a href="https://linuxhandbook.com/date-command/" title="commande date">https://linuxhandbook.com/date-command/</a>

To launch the script, don't forget to allow the execution : <br/> chmod +x ./maintenance.sh

Then go into the folder and execute the script : <br/> ./maintenance.sh

