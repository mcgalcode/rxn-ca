from typing import List

from monty.json import MSONable
import numpy as np
import matplotlib.pyplot as plt

class RecipeStep(MSONable):
    pass

class HeatingStep(RecipeStep):
    """Captures the information about a single step during a heating schedule.
    """

    @classmethod
    def hold(cls, temperature, duration, stage_length = 1):
        return [cls(stage_length, temperature) for _ in range(duration)]
    
    @classmethod
    def sweep(cls, t0, tf, stage_length = 1, temp_step_size = 100):
        if t0 == tf:
            raise ValueError("Initial and final temperatures cannot be the same!")
    
        if tf < t0:
            temp_step_size = -temp_step_size

        temps = [int(s) for s in np.arange(t0,tf,temp_step_size)]
        temps.append(tf)
        return [cls(stage_length, t) for t in temps]

    def __init__(self, duration, temperature):
        self.duration = duration
        self.temperature = temperature

    def as_dict(self):
        d = super().as_dict()
        return {
            **d,
            "duration": self.duration, 
            "temperature": self.temperature
        }
    
class RegrindStep(RecipeStep):
    pass

class HeatingSchedule(MSONable):
    """Captures the information of a heating schedule, e.g. ramping up
    to a particular temperature, holding, and then cooling back down
    """
    
    @classmethod
    def build(cls, *args):
        steps = []
        for step in args:
            if type(step) is list:
                for s in step:
                    steps.append(s)
            else:
                steps.append(step)
        return cls(steps)

    def __init__(self, steps):
        self.steps: List[HeatingStep] = steps

    @property
    def temperature_steps(self):
        return [s for s in self.steps if isinstance(s, HeatingStep)]
    
    @property
    def all_temps(self):
        return list(set([s.temperature for s in self.temperature_steps]))
    
    def temp_at(self, step_idx):
        tallied = 0
        for step in self.temperature_steps:
            tallied += step.duration
            if tallied > step_idx:
                return step.temperature
            
    def get_xy_for_plot(self, step_size: int):
        curr_x = 0
        xs = []
        ys = []


        for step in self.temperature_steps:
            xs.append(curr_x)
            ys.append(step.temperature)
            curr_x += step.duration * step_size
            xs.append(curr_x)
            ys.append(step.temperature)            

        if len(self.temperature_steps) == 1:
            xs.append(self.temperature_steps[0].duration * step_size)
            ys.append(self.temperature_steps[0].temperature)

        return xs, ys        
    
    def plot(self):
        fig, axs = plt.subplots()
        total_length = sum([s.duration for s in self.temperature_steps])

        xs, ys = self.get_xy_for_plot()

        axs.plot(xs, ys)
        axs.hlines(298, -100, total_length + 100, color='r')
        axs.set_title("Heating Schedule")
        axs.set_ylabel("Temperature (K)")
        axs.set_xlabel("Rxn. Coordinate")
    
    def __len__(self):
        return len(self.steps)