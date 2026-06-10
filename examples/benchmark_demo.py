import time
from delineator import delineate, write_outputs, DelineatorConfig


def timer(func):
    """
    Decorator function to time another Python function
    Usage:
       @timer
       def myfunction()

    timer(myfunction)
    """
    def wrapper(*args, **kwargs):
        start_time = time.time()
        result = func(*args, **kwargs)
        end_time = time.time()
        execution_time = round(end_time - start_time, 3)
        print(f"Execution time of '{func.__name__}': {execution_time} seconds")
        return execution_time, result
    return wrapper


@timer
def benchmark_delineate():
    config = DelineatorConfig(
        verbose=True,
        calc_area=False,
        clean=False,
        simplify=False,
        fill=False,
        rivers=False,
        outlets=False,
        data_dir=r"C:\Users\mheberger\Documents\watershed_app\static"
    )
    w, r, o = delineate(17.03333333, 18.67333334, config)


if __name__ == "__main__":
    # For whatever reason, it gets faster after the first call.
    benchmark_delineate()
    benchmark_delineate()
    benchmark_delineate()
