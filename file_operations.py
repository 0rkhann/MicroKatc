"""File helpers for moving generated outputs into result directories"""

import os


def run_directory():
    """Folder all results are written to: the current folder (microkatc.py runs inside results/)"""
    return os.getcwd()


class FileOperations:
    """Moves generated files into output directories under the working directory"""

    @staticmethod
    def move_to_output_directory(output_dir_name, file_name):
        """Moves a file to the specified directory"""
        current_file_path = os.path.join(run_directory(), file_name)
        desired_file_path = os.path.join(run_directory(), output_dir_name, file_name)

        os.makedirs(os.path.dirname(desired_file_path), exist_ok=True)

        if os.path.exists(current_file_path):
            os.rename(current_file_path, desired_file_path)
        else:
            print(f"File not found: {current_file_path}")
