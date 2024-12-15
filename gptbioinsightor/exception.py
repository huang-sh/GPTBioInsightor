from openai import APIStatusError


class ApiKeyMissingError(Exception):
    """Exception raised for missing API_KEY."""
    def __init__(self, message="API_KEY is missing"):
        self.message = message
        super().__init__(self.message)


class ApiBalanceLow(APIStatusError):
    """Exception raised for low balance."""

    def __init__(self, provider=None, status_code=None,message=None):
        if provider == "deepseek" and status_code == 402:
            self.message = "DeepSeek does not have enough quota"
        elif provider == "openai" and status_code == 429:
            self.message = "OpenAI does not have enough quota"
        else:
            #TODO add other provider status code handling
            self.message = message
            # self.message = f"{provider} does not have enough quota"
        super().__init__(message=self.message)
