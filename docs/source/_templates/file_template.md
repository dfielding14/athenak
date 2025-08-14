# File: `[path/to/file.ext]`

## Summary
[1-3 sentence description of what this file contains and its role in the codebase]

## Location
- **Path**: `src/[full/path/to/file.ext]`
- **Module**: [`[ModuleName]`](../modules/[module].md)
- **Include Guard**: `[HEADER_GUARD_NAME_HPP_]` (if applicable)

## Purpose

[Paragraph explaining the file's responsibility, why it exists, and what problem it solves]

## Contents

### Classes/Structures

#### `[ClassName]`
```cpp
class [ClassName] {
  // Key members
};
```
**Purpose**: [What this class represents/does]  
**Inheritance**: `[BaseClass]` (if any)  
**Thread Safety**: `[Thread-safe/Not thread-safe/Device-compatible]`

### Key Functions

#### `[FunctionName]`
```cpp
ReturnType FunctionName(param1, param2);
```
**Purpose**: [What it does]  
**Parameters**:
- `param1`: [Description]
- `param2`: [Description]

**Returns**: [What it returns]  
**Execution Space**: `[Host/Device/Both]`

### Macros/Constants

| Name | Value | Purpose |
|------|-------|---------|
| `[MACRO_NAME]` | `[value]` | [What it's used for] |

### Type Definitions

```cpp
using [TypeAlias] = [ActualType];
```
**Purpose**: [Why this alias exists]

## Dependencies

### Includes
This file includes:
- `[header1.hpp]`: [Why needed]
- `[header2.hpp]`: [Why needed]

### Included By
This file is included by:
- `[file1.cpp]`: [For what purpose]
- `[file2.hpp]`: [For what purpose]

## Device/Host Compatibility

- **Host Functions**: ✅ Available
- **Device Functions**: ✅/❌ [Available/Not available]
- **Kokkos Compatibility**: [Full/Partial/None]

## Usage Examples

### Basic Usage
```cpp
// Example showing typical usage
#include "[thisfile.hpp]"

[ClassName] obj;
obj.Method();
```

### Advanced Usage
```cpp
// More complex example
[code example]
```

## Implementation Notes

### Design Decisions
- [Key design choice and rationale]
- [Another design decision]

### Performance Considerations
- [Memory layout considerations]
- [Computational complexity notes]
- [GPU optimization strategies]

### Algorithms
[Description of key algorithms implemented in this file]

## Potential Issues

⚠️ **Warning**: [Known limitation or gotcha]

⚠️ **Thread Safety**: [Any threading concerns]

⚠️ **Memory**: [Memory management considerations]

## Testing

This file is tested by:
- `[test_file1.cpp]`: [What aspects are tested]
- `[test_file2.cpp]`: [Coverage details]

## Modification History

| Version | Date | Author | Changes |
|---------|------|--------|---------|
| Current | - | - | [Recent changes if known] |

## Related Files

- [`[related1.hpp]`](../reference/[path].md): [Relationship]
- [`[related2.cpp]`](../reference/[path].md): [Relationship]

## See Also

- [Module Documentation](../modules/[module].md)
- [API Reference](../api/[namespace_or_class].rst)
- [Design Pattern](../design/[pattern].md)