# FIR two way crossover network generator

`fir_crossover_taps_generator.py` is intended to be used a quick commnad line crossover network generator. The script generates two text files, `lowpass_crossover.txt` and `highpass_crossover.txt` which are files containing taps of the lowpass and highpass filter respectively. The idea with this is to easily iterate through filter designs, plot them with software of choice and import the taps where needed. The repo also contains `crossover_designer.ipynb` which can be used to see amplitude response of your crossover network.

## Prerequisites
To use the filter generator script, at least [python 3.10](https://www.python.org/downloads/release/python-3100/) is required. Other required packages are:
- [numpy](https://numpy.org/)
- [scipy](https://scipy.org/)

To plot the amplitude spectrum using `crossover_designer.ipynb`, the follwoing packages are required:
- [pandas](https://pandas.pydata.org/)
- [plotly](https://plotly.com/)

On top of any other packages needed for jupyter notebooks to function

## Usage
`fir_crossover_taps_generator.py` is intended to be used as a command line program, and as such comes with a number of flags:

- `-fc`, `--crossover_frequency` - Crossover frequency. If `-fs` is not given, normalization by sampling frequency is assumed. This flag is **mandatory**

- `-Bt`, `--transition_width` - Width of frequency band transition. If `-fs` is not given, normalization by sampling frequency is assumed This flag is **mandatory**

- `-Aa`, `--stopband_attenuation` - Desired stopband attenuation in dB. If `-N` is not given, this flag along width `-Ap` is used to estimate the number of taps needed to reach desired attenuation. This flag is **mandatory if `-N` is not set**. This flag is mandatory when using `-w` option, even when `-N` is set.

- `-Ap`, `--passband-attenuation` - Allowable passband attenuation. If `-N` is not given, this flag along width `-Aa` is used to estimate the number of taps needed to reach desired attenuation. Unlike `-Ap`, this flag is not mandatory unless using `-w`, even when `-N` is set.

- `-N`, `--numtaps` - Number of taps to be used for filter synthesis. This flag is **mandatory if `-Aa` is not set**. This flag overrides the estimation of number of taps if both `-Aa` and it are set. 

- `--numtaps-finder` - Algorithm to be used to determine numtaps. Availale algoritms are [Harris's rule of thumb](https://dsp.stackexchange.com/questions/46303/fred-harris-rule-of-thumb) and [Bellanger's rule of thumb from "Digital Processing of Signals – Theory and Practice"](https://dsp.stackexchange.com/questions/31066/how-many-taps-does-an-fir-filter-need). Defaults to harris if Ap is not give, else bellanger. Is overriden by `-N`. This flag is optional.

- `--fir-algorithm` - Algorithm used to generate taps. Options are firls (`scipy.firls()`) and remez (`scipy.remez()`). Defaults to firls. This flag is optional.

- `-w`, `--weigh-taps` -  Apply error distribution to tap generation algorithm. This option requires that `-Aa` and `-Ap` be set. This flag is optional

- `-o`, `-output-path` - Output tap for filter taps. Defaults to current directory. This flag is optional. 

- `-fs` , `--sampling-frequency` - Sampling frequency. If not specified, inputs are assumed to be normalized by sampling frequency. If specified, input frequencies need to be in Hz. This flag is optional

## Example 
The following example generates a crossover network with a crossover frequency of 270Hz, Transition width of 460Hz, desired attenuation of 80dB, and sampling frequency of 48kHz

```sh
python fir_crossover_taps_generator.py -fc 270 -Bt 460 -Aa 80 -fs 48000
```
The estimated number of taps for this filter is 380. 

## License

[GPL3](https://choosealicense.com/licenses/gpl-3.0/#)