import functools
import inspect

from spatialtis.config import Config


def params(centroid=False, cell_type=False, marker=False, shape=False):
    def use_func(func):
        @functools.wraps(func)
        def wrapper(*args, **kw):
            # to handle default parameters
            sig = inspect.signature(func)
            func_params = sig.parameters.keys()
            user_kw_keys = kw.keys()

            if "exp_obs" in user_kw_keys:
                exp_obs = kw["exp_obs"]
                if isinstance(exp_obs, (str, int, float)):
                    kw["exp_obs"] = [exp_obs]
                else:
                    kw["exp_obs"] = list(exp_obs)

            def set_kwargs(kw_name, default_value, error=True):
                if kw_name in func_params:
                    if kw_name not in user_kw_keys:  # if not in user input, we get the default value
                        if default_value is not None:
                            kw[kw_name] = default_value
                        else:
                            if error:
                                raise ValueError(f"Please set Config.{kw_name} or pass `{kw_name}=`")

            set_kwargs("exp_obs", Config.exp_obs)
            set_kwargs("mp", Config.mp, error=False)
            set_kwargs("cell_type_key", Config.cell_type_key) if cell_type else None
            set_kwargs("centroid_key", Config.centroid_key) if centroid else None
            set_kwargs("shape_key", Config.shape_key) if shape else None
            set_kwargs("marker_key", Config.marker_key) if marker else None

            return func(*args, **kw)
        return wrapper
    return use_func
