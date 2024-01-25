# README.md

## _0_main.py

When you run `_0_main.py`, it will display the growth series and other data of SACADA #1.
To calculate for a different example, please modify the `data` around line 75 and `s` around line 98.

The folder `SACADApy` contains several data of SACADA (please refer to [HKGP16] and the website [sacada.info](https://www.sacada.info)). We have transformed original data mainly to reduce the kinds of vertices.

The folder `tilingdata` contains other several examples. To calculate for your own examples, follow the instructions below.

1. Create input data within the folder `tilingdata`. For the input format, please refer [the main page](../README.md).
1. Import the created input data at the beginning of `_0_main.py` and modify the variables `data` (and `s` if needed).

## reference

[1] T. Inoue and Y. Nakamura, Ehrhart theory on periodic graphs, [arXiv: 2305.08177v1](https://arxiv.org/abs/2305.08177)

[HKGP16] R. Hoffmann, A. A. Kabanov, A. A. Golov, and D. M. Proserpio, Homo Citans and Carbon Allotopes: For an Ethics of Citation, Angewandte Chemie International Edition 55 (2016), no. 37, 10962–10976.

## License

See the [LICENSE](LICENSE.md) file for license rights and limitations (MIT).
