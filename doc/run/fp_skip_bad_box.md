# Configuration Filtering with fp_skip_bad_box

The `fp_skip_bad_box` parameter provides an efficient way to filter out configurations with problematic simulation box geometries before performing expensive first-principles calculations. This can significantly improve the efficiency of DP-GEN workflows by avoiding calculations on physically unreasonable structures.

## Purpose

During molecular dynamics simulations used for model deviation calculations, some configurations may develop problematic box geometries that can cause issues in first-principles calculations:

- Extremely elongated or compressed cells that violate periodic boundary conditions
- Highly tilted triclinic cells that may be difficult to converge
- Poor aspect ratios that lead to numerical instabilities
- Excessive wrapping of the simulation cell

The `fp_skip_bad_box` feature automatically identifies and skips such configurations, improving both efficiency and robustness of the DP-GEN workflow.

## Parameter Format

The parameter accepts a semicolon-separated string of colon-separated key-value pairs:

```json
"fp_skip_bad_box": "criterion1:threshold1;criterion2:threshold2;..."
```

## Available Criteria

### 1. length_ratio

**Purpose**: Controls the maximum ratio of cell edge lengths (longest/shortest).

**Formula**: `max(|a|, |b|, |c|) / min(|a|, |b|, |c|)`

**Description**: Prevents extremely elongated or compressed cells by ensuring the cell dimensions are reasonably balanced.

**Suggested values**: 3-5

**Example**: `"length_ratio:3"` will skip configurations where one cell edge is more than 3 times longer than the shortest edge.

### 2. height_ratio

**Purpose**: Controls the ratio of the maximum cell edge length to the minimum face-to-face distance.

**Formula**: `max(|a|, |b|, |c|) / min(face_distances)`

**Description**: Prevents "needle-like" or "pancake-like" cells where the cell is very thin in one direction relative to its longest dimension.

**Suggested values**: 3-5

**Example**: `"height_ratio:3"` ensures reasonable thickness in all directions relative to the cell size.

### 3. wrap_ratio

**Purpose**: Controls triclinic wrapping by limiting off-diagonal cell matrix elements.

**Formula**: `max(|b_x/a_x|, |c_y/b_y|, |c_x/a_x|)`

**Description**: Prevents excessive shearing and wrapping in triclinic cells that can cause visualization and convergence issues.

**Suggested values**: 0.3-0.7

**Example**: `"wrap_ratio:0.5"` limits how much the cell can be sheared from orthogonal.

### 4. tilt_ratio

**Purpose**: Controls cell tilt angles by limiting the ratio of off-diagonal to diagonal elements.

**Formula**: `max(|b_x/b_y|, |c_y/c_z|, |c_x/c_z|)`

**Description**: Prevents extreme tilt angles that can cause issues with periodic boundary conditions and force calculations.

**Suggested values**: 0.3-0.7

**Example**: `"tilt_ratio:0.5"` ensures moderate tilt angles in triclinic cells.

## Usage Examples

### Basic Example
```json
{
  "fp_skip_bad_box": "length_ratio:3;height_ratio:3"
}
```

### Comprehensive Filtering
```json
{
  "fp_skip_bad_box": "length_ratio:3;height_ratio:3;wrap_ratio:0.5;tilt_ratio:0.5"
}
```

### Conservative Filtering (allows more deformed cells)
```json
{
  "fp_skip_bad_box": "length_ratio:5;height_ratio:5;wrap_ratio:0.7;tilt_ratio:0.7"
}
```

### Strict Filtering (rejects more deformed cells)
```json
{
  "fp_skip_bad_box": "length_ratio:2;height_ratio:2;wrap_ratio:0.3;tilt_ratio:0.3"
}
```

## Best Practices

1. **Start conservative**: Begin with lenient thresholds and tighten them if needed based on your system's behavior.

2. **System-dependent tuning**: Different materials and conditions may require different thresholds:
   - Fluid systems: More lenient wrapping and tilt ratios
   - Crystal systems: Stricter length and height ratios
   - High-pressure simulations: May need adjusted thresholds

3. **Monitor filtering statistics**: Check DP-GEN logs to see how many configurations are being filtered and adjust accordingly.

4. **Balance efficiency vs completeness**: Overly strict filtering may miss important rare configurations, while too lenient filtering may waste computational resources.

## Technical Details

The filtering is applied in the `make_fp` stage of DP-GEN, after model deviation calculations but before first-principles calculations. This placement ensures:

- Maximum computational savings (FP calculations are the most expensive)
- Preservation of exploration diversity (MD simulations run to completion)
- Accurate model deviation statistics (all MD configurations are evaluated)

The cell matrix format follows the standard convention where the cell vectors are rows of a 3×3 matrix:
```
cell = [[a_x, a_y, a_z],
        [b_x, b_y, b_z], 
        [c_x, c_y, c_z]]
```

## Common Issues and Solutions

### Too many configurations filtered
- Increase threshold values
- Check if your MD simulations are producing reasonable geometries
- Consider adjusting MD parameters (temperature, pressure, timestep)

### First-principles calculations failing
- Decrease threshold values to filter more aggressively
- Add criteria that weren't previously used
- Check specific failure modes to identify appropriate criteria

### Unexpected filtering behavior
- Verify the cell matrix convention matches your system
- Check that the format string is correctly parsed
- Review DP-GEN logs for detailed filtering statistics