import uuid


def output_file_name(*,
                     extension='.res',
                     file_name_length: int = 5) -> str:
    uuid_ = str(uuid.uuid4())
    base_name = uuid_[:file_name_length]
    return ''.join([base_name, extension])
