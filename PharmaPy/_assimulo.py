"""Lazy internal constructors for the optional Assimulo solver backend."""

from importlib import import_module
from typing import Any

_INSTALL_MESSAGE = (
    "Assimulo is required for solver-backed PharmaPy simulations. "
    "Install Assimulo 3.4.3 from conda-forge using environment.yml or "
    "provide a compatible local installation."
)


def _load_assimulo_symbol(module_name: str, symbol_name: str) -> Any:
    """Load a symbol from the optional Assimulo dependency.

    Parameters
    ----------
    module_name : str
        Fully qualified Assimulo module name.
    symbol_name : str
        Attribute to load from the module.

    Returns
    -------
    Any
        Requested Assimulo class or callable.

    Raises
    ------
    ImportError
        If Assimulo or one of its required runtime components cannot be
        imported.
    """
    try:
        module = import_module(module_name)
    except ImportError as exc:
        raise ImportError(_INSTALL_MESSAGE) from exc

    try:
        return getattr(module, symbol_name)
    except AttributeError as exc:
        raise ImportError(
            f"{_INSTALL_MESSAGE} The installed version does not provide "
            f"{module_name}.{symbol_name}."
        ) from exc


def _construct_assimulo_object(
    module_name: str,
    symbol_name: str,
    *args: Any,
    **kwargs: Any,
) -> Any:
    """Construct an object from the optional Assimulo dependency.

    Parameters
    ----------
    module_name : str
        Fully qualified Assimulo module name.
    symbol_name : str
        Constructor name within ``module_name``.
    *args : Any
        Positional arguments forwarded unchanged to the constructor.
    **kwargs : Any
        Keyword arguments forwarded unchanged to the constructor.

    Returns
    -------
    Any
        Constructed Assimulo problem or solver object.

    Raises
    ------
    ImportError
        If the requested Assimulo constructor is unavailable.
    """
    constructor = _load_assimulo_symbol(module_name, symbol_name)
    return constructor(*args, **kwargs)


def CVode(*args: Any, **kwargs: Any) -> Any:
    """Construct an Assimulo CVode solver on first use.

    Parameters
    ----------
    *args : Any
        Positional arguments forwarded to ``assimulo.solvers.CVode``.
    **kwargs : Any
        Keyword arguments forwarded to ``assimulo.solvers.CVode``.

    Returns
    -------
    Any
        Assimulo CVode solver instance.

    Raises
    ------
    ImportError
        If the optional Assimulo backend is unavailable.
    """
    return _construct_assimulo_object("assimulo.solvers", "CVode", *args, **kwargs)


def IDA(*args: Any, **kwargs: Any) -> Any:
    """Construct an Assimulo IDA solver on first use.

    Parameters
    ----------
    *args : Any
        Positional arguments forwarded to ``assimulo.solvers.IDA``.
    **kwargs : Any
        Keyword arguments forwarded to ``assimulo.solvers.IDA``.

    Returns
    -------
    Any
        Assimulo IDA solver instance.

    Raises
    ------
    ImportError
        If the optional Assimulo backend is unavailable.
    """
    return _construct_assimulo_object("assimulo.solvers", "IDA", *args, **kwargs)


def Explicit_Problem(*args: Any, **kwargs: Any) -> Any:
    """Construct an Assimulo explicit problem on first use.

    Parameters
    ----------
    *args : Any
        Positional arguments forwarded to
        ``assimulo.problem.Explicit_Problem``.
    **kwargs : Any
        Keyword arguments forwarded to
        ``assimulo.problem.Explicit_Problem``.

    Returns
    -------
    Any
        Assimulo explicit-problem instance.

    Raises
    ------
    ImportError
        If the optional Assimulo backend is unavailable.
    """
    return _construct_assimulo_object(
        "assimulo.problem", "Explicit_Problem", *args, **kwargs
    )


def Implicit_Problem(*args: Any, **kwargs: Any) -> Any:
    """Construct an Assimulo implicit problem on first use.

    Parameters
    ----------
    *args : Any
        Positional arguments forwarded to
        ``assimulo.problem.Implicit_Problem``.
    **kwargs : Any
        Keyword arguments forwarded to
        ``assimulo.problem.Implicit_Problem``.

    Returns
    -------
    Any
        Assimulo implicit-problem instance.

    Raises
    ------
    ImportError
        If the optional Assimulo backend is unavailable.
    """
    return _construct_assimulo_object(
        "assimulo.problem", "Implicit_Problem", *args, **kwargs
    )
