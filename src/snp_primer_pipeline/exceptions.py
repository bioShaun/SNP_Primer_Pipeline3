"""Custom exceptions for SNP primer pipeline."""


class PipelineError(Exception):
    """Base exception for all pipeline errors."""
    pass


class ParseError(PipelineError):
    """Exception raised during input parsing."""
    
    def __init__(self, message: str, line_number: int = None, line_content: str = None):
        self.line_number = line_number
        self.line_content = line_content
        
        if line_number is not None:
            message = f"Line {line_number}: {message}"
        if line_content is not None:
            message = f"{message} (content: {line_content[:50]}...)"
            
        super().__init__(message)


class BlastError(PipelineError):
    """Exception raised during BLAST execution or parsing."""
    
    def __init__(self, message: str, command: str = None, return_code: int = None):
        self.command = command
        self.return_code = return_code
        
        if command is not None:
            message = f"BLAST command failed: {command}\n{message}"
        if return_code is not None:
            message = f"{message} (exit code: {return_code})"
            
        super().__init__(message)


class AlignmentError(PipelineError):
    """Exception raised during multiple sequence alignment."""
    
    def __init__(self, message: str, snp_name: str = None, sequence_count: int = None):
        self.snp_name = snp_name
        self.sequence_count = sequence_count
        
        if snp_name is not None:
            message = f"Alignment failed for SNP {snp_name}: {message}"
        if sequence_count is not None:
            message = f"{message} ({sequence_count} sequences)"
            
        super().__init__(message)


class PrimerDesignError(PipelineError):
    """Exception raised during primer design."""
    
    def __init__(self, message: str, snp_name: str = None, design_type: str = None):
        self.snp_name = snp_name
        self.design_type = design_type
        
        if snp_name is not None:
            message = f"Primer design failed for SNP {snp_name}: {message}"
        if design_type is not None:
            message = f"{design_type} {message}"
            
        super().__init__(message)


class ConfigurationError(PipelineError):
    """Exception raised for configuration errors."""
    
    def __init__(self, message: str, config_file: str = None, parameter: str = None):
        self.config_file = config_file
        self.parameter = parameter
        
        if config_file is not None:
            message = f"Configuration error in {config_file}: {message}"
        if parameter is not None:
            message = f"{message} (parameter: {parameter})"
            
        super().__init__(message)