from os import PathLike
from typing import TYPE_CHECKING, Protocol, Union, runtime_checkable

# Some classes are declared as generic in stubs, but not at runtime.
# In Python 3.9 and earlier, os.PathLike is not subscriptable, results in runtime error
if TYPE_CHECKING:
    from builtins import ellipsis

    Index = ellipsis | slice | int
    PathLikeStr = PathLike[str]
else:
    PathLikeStr = PathLike


@runtime_checkable
class GeneralPath(Protocol):
    """A protocol to handle paths that can be either local or S3 paths."""

    def parent(self): ...

    def suffix(self): ...

    def resolve(self): ...

    def exists(self): ...

    def read_text(self): ...

    def __truediv__(self, other): ...

    def __str__(self) -> str: ...

    def __fspath__(self) -> str:
        return str(self)


PathOrStr = Union[str, PathLikeStr, GeneralPath]
