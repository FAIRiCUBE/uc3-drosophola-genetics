from measurer import Measurer
from types import ModuleType


def your_function():
    # Path where the data are stored (the use of the disk in this path is measured).
    # Use '/' to measure the entire disk.
    data_path = '/'
    measurer = Measurer()
    tracker = measurer.start(data_path=data_path)
    # example -> shape = [5490, 2170]
    shape = []

    '''
        # write here your code...
    '''

    # it is very important to use program_path = __file__
    measurer.end(tracker=tracker,
                 shape=shape,
                 libraries=[v.__name__ for k, v in globals().items() if type(v) is ModuleType and not k.startswith('__')],
                 data_path=data_path,
                 program_path=__file__,
                 variables=locals(),
                 csv_file='benchmarks.csv')
