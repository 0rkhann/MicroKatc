import os

CURRENT_DIRECTORY = os.getcwd()

class FileOperations:
    @staticmethod
    def move_to_output_directory(output_dir_name, file_name):
        """Moves a file to the specified directory"""
        current_file_path = os.path.join(CURRENT_DIRECTORY, file_name)
        desired_file_path = os.path.join(CURRENT_DIRECTORY, output_dir_name, file_name)

        os.makedirs(os.path.dirname(desired_file_path), exist_ok = True)

        if os.path.exists(current_file_path):
            os.rename(current_file_path, desired_file_path)
        else:
            print(f"File not found: {current_file_path}")