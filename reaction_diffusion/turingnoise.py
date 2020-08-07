

import time

import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
from scipy.signal import convolve2d
from pathlib import Path
import imageio

class DiffMatrix:
    def __init__(self, size=300, start_template="square", reaction = 'Coral'):
        
        self.start_template = start_template
        self.size = size
        self.reaction = reaction

        self.start_time = int(time.time())
        Path(f'matrix_plots/{self.start_time}').mkdir()

        self.matrix_shape = (int(self.size), int(self.size))
        self.matrix_a = self.initiate_matrix()
        self.matrix_b = self.initiate_matrix()

        ### Reaction caracteristics ###
        if self.reaction == "Coral":
            self.diff_a = 1
            self.diff_b = 0.5
            self.feed_rate = 0.055
            self.kill_rate = 0.062
            self.timestep = 1
            self.laplace_window = 3

        self.prepare_matrix()
        self.initialize_sum_matricies()

    def initiate_matrix(self):
        matrix = np.zeros(self.matrix_shape)
        return matrix

    def prepare_matrix(self):
        self.matrix_a = np.ones(self.matrix_shape)
        self.matrix_b = np.zeros(self.matrix_shape)

        if self.start_template == "line":
            start_b = int(0.4 * self.size)
            end_b = int(0.6 * self.size)
            self.matrix_b[start_b:end_b] = 1

        elif self.start_template == "square":
            start_b = int(0.45 * self.size)
            end_b = int(0.55 * self.size)
            self.matrix_b[start_b:end_b, start_b:end_b] = 1

    def initialize_sum_matricies(self,):
        self.matrix_a_feed = self.initiate_matrix()
        self.matrix_a_react = self.initiate_matrix()
        self.matrix_a_diffuse = self.initiate_matrix()

        self.matrix_b_feed = self.initiate_matrix()
        self.matrix_b_react = self.initiate_matrix()
        self.matrix_b_diffuse = self.initiate_matrix()

    def feed(self):
        # create a up until 1
        self.matrix_a_feed = self.feed_rate * (1 - self.matrix_a)
        # bestroy b if any
        self.matrix_b_feed = (self.feed_rate + self.kill_rate) * self.matrix_b

    def react(self):
        # consume a and turn it into b
        
        reaction_rate = self.matrix_a * self.matrix_b * self.matrix_b
        
        self.matrix_a_react = reaction_rate
        self.matrix_b_react = reaction_rate

    def diffuse(self):
        diffusion_kernel = np.array([
            [0.05, 0.2, 0.05], 
            [0.2, -1, 0.2], 
            [0.05, 0.2, 0.05]])

        self.matrix_a_diffuse = self.diff_a * convolve2d(self.matrix_a, diffusion_kernel, boundary="wrap",mode='same')
        self.matrix_b_diffuse = self.diff_b * convolve2d(self.matrix_b, diffusion_kernel, boundary="wrap",mode='same')

    def plot_matrix_to_file(self,i):
        fig,axes = plt.subplots(1,2)
        axes[0].imshow(dm.matrix_a, cmap='binary', interpolation='nearest',vmin=0, vmax=1)
        axes[0].set_title(f'Activator at {i}')
        axes[1].set_title(f'Blocker at {i}')
        axes[1].imshow(dm.matrix_b, cmap='binary', interpolation='nearest',vmin=0, vmax=1)
        plt.savefig(f'matrix_plots/{self.start_time}/matrix_{i}.png')


    def print_matrix(self):
        print(self.matrix_a)
        print('====================')
        print(self.matrix_b)

    def print_matrix_change(self):
        print('\t### matrix_a_diffuse')
        print(self.matrix_a_diffuse)
        print('\t### matrix_a_react')
        print(self.matrix_a_react)
        print('\t### matrix_a_feed')
        print(self.matrix_a_feed)
        print('====================')
        print('\t### matrix b diffuse')
        print(self.matrix_b_diffuse)
        print('\t### matrix_b_react')
        print(self.matrix_b_react)
        print('\t### matrix_b_feed')
        print(self.matrix_b_feed)
        

    def _next(self):
        # print(f"Mean a: {dm.matrix_a.mean()}Mean b:{dm.matrix_b.mean()}")
        # print(f"Max a: {dm.matrix_b.max()}Max b:{dm.matrix_b.max()}")

        self.diffuse()
        self.feed()
        self.react()

        self.matrix_a = self.matrix_a + (self.matrix_a_diffuse - self.matrix_a_react + self.matrix_a_feed)
        self.matrix_b = self.matrix_b + (self.matrix_b_diffuse + self.matrix_b_react - self.matrix_b_feed)
        
        self.initialize_sum_matricies()


if __name__ == "__main__":
    print_every_n = 100
    dm = DiffMatrix(400, start_template='square',reaction='Mitosis')

    for i in tqdm(range(20000)):
        
        dm._next()
        if i % print_every_n == 0:
            dm.plot_matrix_to_file(i)
    