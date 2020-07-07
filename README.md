# Ion optics

### Optimization of Twiss parameters to influence the emissivity of an ion beam.

<p>Programs developed during an internship in ion optics to discover the notion of emissivity, Twiss parameters and "sigma" matrices for ion beams passing through different optical systems (drift, thin lens, Einzel lens, magnetic dipole, electrostatic quadrupole).</p>

<p>The parameters are given in "parameters.py", which allows to change the values of the optical elements. The Twiss parameter "gamma" is calculated in "gamma.py".</p>

<p>The "trash cleanup" part allows to remove (or not) files that have been in the trash for more than 1 week. But you can change this value with the corresponding macro. Explanations here: <a href="https://linuxhandbook.com/date-command/" title="commande date">https://linuxhandbook.com/date-command/</a>

To launch the script, don't forget to allow the execution : <br/> chmod +x ./maintenance.sh

Then go into the folder and execute the script : <br/> ./maintenance.sh

And voilà, everything is done by itself:
<b><ol>
    <li>Updating</li>
    <li>Autoremove/clean</li>
    <li>"Localepurge" junk files according to your language</li>
    <li>Purging deleted packages</li>
    <li>Purge useless kernels and put them in "manuals".</li>
    <li>Snapshot update</li>
    <li>Resolution of dependencies not present</li>
    <li>Cleaning the basket</li>
</ol>
</b>