# Terrain amplification (or super-resolution)
Terrain amplification implementation in Matlab as described in:
Sparse representation of terrains for procedural modeling
Eric Guérin, Julie Digne, Eric Galin, Adrien Peytavie
Computer Graphics Forum. Volume 35, number 2. Proceedings of Eurographics 2016, Lisbon.

## Dependencies
You must install ompbox and ksvdbox in order to execute this code.

## Usage
```Matlab
terrain_super_resolution( factor, exemplarhr, inputterrain, output, masksize, offset_analysis, offset_synthesis)
```

Where:
- `factor` is the amplification factor
- `exemplarhr` : the high resolution exemplar filename
- `inputterrain` : the input terrain filename
- `output` : the output terrain filename
- `masksize` : the diameter of the patch on the low resolution input terrain and exemplar
- `offset_analysis` : offset of the patch extracted from the low resolution exemplar
- `offset_synthesis` : offset that will be used for the synthesis

## Example
```Matlab
terrain_super_resolution(4,'grandcanyonhr.png','sketchlr.png',16,8,8);
```

## Warning
When  using too small values for offset or mask size, or too large images as exemplar, you could have memory issues.

## Reference
If you use this code, you must reference this article:

```
@article {guerin2016,
	author = {Guérin, Eric and Digne, Julie and Galin, 
		Eric and Peytavie, Adrien},
	title = {Sparse representation of terrains for 
		procedural modeling},
	journal = {Computer Graphics Forum 
		(proceedings of Eurographics 2016)},
	volume = {35},
	number = {2},
	issn = {1467-8659},
	url = {http://dx.doi.org/10.1111/cgf.12821},
	doi = {10.1111/cgf.12821},
	pages = {177--187},
	year = {2016},
}
```
