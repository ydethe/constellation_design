from pdm.backend.hooks.base import Context
from pybind11.setup_helpers import Pybind11Extension


def pdm_build_update_setup_kwargs(context: Context, setup_kwargs: dict):
    if context.target != "sdist":
        ext_modules = [
            Pybind11Extension(
                name="constellation_design.simulate",
                sources=[
                    "src/optim_fn.cpp",
                    "src/satellite.cpp",
                    "src/simulate.cpp",
                    "src/utils.cpp",
                ],
                include_dirs=["/usr/local/include/eigen3", "src"],
                extra_compile_args=["-Os", "-std=c++11"],
                extra_link_args=["-lm", "-lnlopt", "-lgsl"],
            )
        ]
        setup_kwargs.update(
            setup_requires=["pybind11"],
            ext_modules=ext_modules,
        )
