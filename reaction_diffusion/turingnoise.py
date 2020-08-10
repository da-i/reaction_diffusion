import time

import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
from scipy.signal import convolve2d
from pathlib import Path
import imageio


class DiffMatrix:
    def __init__(self, size=300, start_template="square", reaction="Coral"):

        self.start_template = start_template
        self.size = size
        self.reaction = reaction

        self.start_time = int(time.time())
        self.img_dir = Path(f"matrix_plots/{self.start_time}")
        self.img_dir.mkdir()

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

        elif self.reaction == "Mitosis":
            self.diff_a = 1
            self.diff_b = 0.5
            self.feed_rate = 0.0367
            self.kill_rate = 0.0649
            self.timestep = 1
            self.laplace_window = 3

        elif self.reaction == "test":
            self.diff_a = 1
            self.diff_b = 0.5
            # y axis scaled
            self.feed_rate = self.make_range_array(0.001, 0.1, self.size)
            self.feed_rate = np.transpose(self.feed_rate)
            #x axis scaled
            self.kill_rate = self.make_range_array(0.05, 0.065, self.size)

            self.timestep = 1
            self.laplace_window = 3

        self.prepare_matrix()
        self.initialize_sum_matricies()

    @staticmethod
    def make_range_array(start, stop, length):
        step = (stop - start) / length
        array1d = np.arange(start, stop, step)[:length]
        # turn the 1d array into a 2d array by repeating it n times
        return np.tile(array1d, (length, 1))

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

        elif self.start_template == "smallsquare":
            start_b = int(0.48 * self.size)
            end_b = int(0.52 * self.size)
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
        diffusion_kernel = np.array(
            [[0.05, 0.2, 0.05], [0.2, -1, 0.2], [0.05, 0.2, 0.05]]
        )

        self.matrix_a_diffuse = self.diff_a * convolve2d(
            self.matrix_a, diffusion_kernel, boundary="wrap", mode="same"
        )
        self.matrix_b_diffuse = self.diff_b * convolve2d(
            self.matrix_b, diffusion_kernel, boundary="wrap", mode="same"
        )

    def plot_matrix_to_file(self, i, style="single"):
        i_fmt = "{:09d}".format(i)

        if style == "single":
            fig = plt.imshow(dm.matrix_a, cmap="binary", interpolation="spline16")
            plt.title(f"Activator at {i}")
            plt.axis("off")

        elif style == "double":
            fig, axes = plt.subplots(1, 2)
            axes[0].imshow(
                dm.matrix_a, cmap="binary", interpolation="nearest", vmin=0, vmax=1
            )
            axes[0].set_title(f"Activator at {i}")
            axes[1].set_title(f"Blocker at {i}")
            axes[1].imshow(
                dm.matrix_b, cmap="binary", interpolation="nearest", vmin=0, vmax=1
            )

        plt.savefig(f"matrix_plots/{self.start_time}/{i_fmt}.png")
        plt.close()

    def print_matrix(self):
        print(self.matrix_a)
        print("====================")
        print(self.matrix_b)

    def print_matrix_change(self):
        print("\t### matrix_a_diffuse")
        print(self.matrix_a_diffuse)
        print("\t### matrix_a_react")
        print(self.matrix_a_react)
        print("\t### matrix_a_feed")
        print(self.matrix_a_feed)
        print("====================")
        print("\t### matrix b diffuse")
        print(self.matrix_b_diffuse)
        print("\t### matrix_b_react")
        print(self.matrix_b_react)
        print("\t### matrix_b_feed")
        print(self.matrix_b_feed)

    def generate_gif(self):
        filenames = self.img_dir.glob("*.png")
        filenames = [x for x in filenames]
        filenames.sort()
        with imageio.get_writer(
            f"{self.img_dir}/{self.reaction}.gif", mode="I"
        ) as writer:
            for filename in filenames:
                image = imageio.imread(filename)
                writer.append_data(image)

    def _next(self):

        # calculate changes to system
        self.diffuse()
        self.feed()
        self.react()

        # Update matricies
        self.matrix_a += (
            self.matrix_a_diffuse - self.matrix_a_react + self.matrix_a_feed
        )
        self.matrix_b += (
            self.matrix_b_diffuse + self.matrix_b_react - self.matrix_b_feed
        )

        self.initialize_sum_matricies()


if __name__ == "__main__":
    print_every_n = 50
    dm = DiffMatrix(400, start_template="square", reaction="test")
    # dm.fill_matrix_random()

    for i in tqdm(range(20000)):

        dm._next()
        if i % print_every_n == 0:
            dm.plot_matrix_to_file(i)

    dm.generate_gif()
