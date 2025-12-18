def bounds_contains(
    outer_bounds: tuple[float, float, float, float] | dict[str, float],
    inner_bounds: tuple[float, float, float, float] | dict[str, float],
    buffer: float = 0.0,
) -> bool:
    """Check if outer bounds completely contain inner bounds.

    Parameters
    ----------
    outer_bounds : tuple or dict
        Outer bounds as (west, south, east, north) or dict with those keys.
    inner_bounds : tuple or dict
        Inner bounds as (west, south, east, north) or dict with those keys.
    buffer : float, optional
        Buffer distance to require around inner bounds (in same units).
        Default is 0.0 (exact containment).

    Returns
    -------
    bool
        True if outer bounds completely contain inner bounds (with buffer).

    Examples
    --------
    >>> # Check if UNR grid covers DISP frame
    >>> unr_bounds = (-97, 28, -93, 32)
    >>> disp_bounds = (-96, 29, -94, 31)
    >>> bounds_contains(unr_bounds, disp_bounds)
    True

    >>> # With buffer requirement
    >>> bounds_contains(unr_bounds, disp_bounds, buffer=0.5)
    True

    >>> # Dict format
    >>> unr_bounds = {"west": -97, "south": 28, "east": -93, "north": 32}
    >>> disp_bounds = {"west": -96, "south": 29, "east": -94, "north": 31}
    >>> bounds_contains(unr_bounds, disp_bounds)
    True

    >>> # Does not contain
    >>> small_bounds = (-95, 29, -94, 30)
    >>> large_bounds = (-96, 28, -93, 31)
    >>> bounds_contains(small_bounds, large_bounds)
    False

    """
    # Parse outer bounds
    if isinstance(outer_bounds, dict):
        o_west = outer_bounds["west"]
        o_south = outer_bounds["south"]
        o_east = outer_bounds["east"]
        o_north = outer_bounds["north"]
    else:
        o_west, o_south, o_east, o_north = outer_bounds

    # Parse inner bounds
    if isinstance(inner_bounds, dict):
        i_west = inner_bounds["west"]
        i_south = inner_bounds["south"]
        i_east = inner_bounds["east"]
        i_north = inner_bounds["north"]
    else:
        i_west, i_south, i_east, i_north = inner_bounds

    # Check containment with buffer
    return (
        o_west <= i_west - buffer
        and o_south <= i_south - buffer
        and o_east >= i_east + buffer
        and o_north >= i_north + buffer
    )


def check_bounds_coverage(
    outer_bounds: tuple[float, float, float, float] | dict[str, float],
    inner_bounds: tuple[float, float, float, float] | dict[str, float],
    buffer: float = 0.0,
) -> dict[str, bool | dict[str, float]]:
    """Check bounds coverage with detailed gap information.

    Parameters
    ----------
    outer_bounds : tuple or dict
        Outer bounds as (west, south, east, north) or dict.
    inner_bounds : tuple or dict
        Inner bounds as (west, south, east, north) or dict.
    buffer : float, optional
        Buffer distance required. Default is 0.0.

    Returns
    -------
    dict
        Dictionary with:
        - "contains": bool, True if outer contains inner (with buffer)
        - "gaps": dict, Gap distances by direction (negative = covered)

    Examples
    --------
    >>> unr_bounds = (-97, 28, -93, 32)
    >>> disp_bounds = (-96, 29, -94, 31)
    >>> result = check_bounds_coverage(unr_bounds, disp_bounds, buffer=0.5)
    >>> result["contains"]
    True
    >>> result["gaps"]
    {'west': -0.5, 'south': -0.5, 'east': -0.5, 'north': -0.5}

    >>> # Insufficient coverage
    >>> small_bounds = (-95.5, 29, -94, 30)
    >>> result = check_bounds_coverage(small_bounds, disp_bounds)
    >>> result["contains"]
    False
    >>> result["gaps"]
    {'west': 0.5, 'south': 0.0, 'east': 0.0, 'north': 1.0}

    """
    # Parse bounds
    if isinstance(outer_bounds, dict):
        o_west = outer_bounds["west"]
        o_south = outer_bounds["south"]
        o_east = outer_bounds["east"]
        o_north = outer_bounds["north"]
    else:
        o_west, o_south, o_east, o_north = outer_bounds

    if isinstance(inner_bounds, dict):
        i_west = inner_bounds["west"]
        i_south = inner_bounds["south"]
        i_east = inner_bounds["east"]
        i_north = inner_bounds["north"]
    else:
        i_west, i_south, i_east, i_north = inner_bounds

    # Calculate gaps (positive = missing coverage, negative = extra coverage)
    gaps = {
        "west": (i_west - buffer) - o_west,
        "south": (i_south - buffer) - o_south,
        "east": o_east - (i_east + buffer),
        "north": o_north - (i_north + buffer),
    }

    # Check if all gaps are negative or zero (fully contained)
    contains = all(gap <= 0 for gap in gaps.values())

    return {
        "contains": contains,
        "gaps": gaps,
    }
