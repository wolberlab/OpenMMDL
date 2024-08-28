import re

class ConfigParser:
    def __init__(self, config_file_path):
        self.config_file_path = config_file_path
        self.variables = {}
        self._parse_config_file()
        self._set_attributes()

    def _parse_config_file(self):
        # Define a regular expression to match variable assignments
        variable_pattern = re.compile(r'^\s*(\w+)\s*=\s*(.+?)\s*$')

        # Read the file and extract the variables
        with open(self.config_file_path, 'r') as file:
            for line in file:
                match = variable_pattern.match(line)
                if match:
                    var_name = match.group(1)
                    var_value = match.group(2).strip('\'"')  # Remove surrounding quotes if any
                    # Store the variable and its value in the dictionary
                    self.variables[var_name] = var_value

    def _set_attributes(self):
        for var_name, var_value in self.variables.items():
            setattr(self, var_name, var_value)

    def get_variable(self, var_name, default=None):
        return self.variables.get(var_name, default)
    
    def __getattr__(self, name):
        # This method is called when the accessed attribute is not found
        return self.variables.get(name, None)

    def print_variables(self):
        for var_name, var_value in self.variables.items():
            print(f'{var_name} = {var_value}')