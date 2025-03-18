from measurer import Measurer
from types import ModuleType
import boto3


aws_session = boto3.Session(
    aws_access_key_id='<your_access_key_id>',
    aws_secret_access_key='<your_secret_access_key>'
)


def your_function():
    # Use the name of the AWS S3 bucket
    data_path = 'AWS-S3-bucket-name'
    measurer = Measurer()
    tracker = measurer.start(data_path=data_path, aws_s3=True, aws_session=aws_session)
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
                 csv_file='benchmarks.csv',
                 aws_s3=True,
                 aws_session=aws_session)
