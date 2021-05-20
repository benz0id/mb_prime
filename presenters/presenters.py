from abc import ABC, abstractmethod


class Presenter(ABC):
    """An abstract class outlining the basic functionality of any presenter."""

    def __init__(self):
        pass

    @abstractmethod
    def print(self, txt: str) -> None:
        """Displays the text in some manner"""
        pass

class ConsolerPresenter(Presenter):
    """A simple presenter designed to display text to the console."""

    # override
    def print(self, txt: str) -> None:
        """Displays text to the console."""
        print(txt)
        return None

