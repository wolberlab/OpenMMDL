import os
import pylab
from PIL import Image


class FigureMerger:
    """
    Handles the creation and merging of binding mode figures with corresponding legends.

    Attributes
    ----------
    binding_mode : str
        Name of the binding mode for which figures are created.
    occurrence_percent : float
        Occurrence percentage of the binding mode.
    split_data : list of str
        Interaction descriptors used for generating the figure legend.
    merged_image_paths : list of str
        List storing paths to the output merged images.
    """
    def __init__(
        self, binding_mode, occurrence_percent, split_data, merged_image_paths
    ):
        self.binding_mode = binding_mode
        self.occurrence_percent = occurrence_percent
        self.split_data = split_data
        self.merged_image_paths = merged_image_paths

    def create_and_merge_images(self):
        """
        Create and merge images to generate a legend for binding modes.

        Returns
        -------
        list of str
            Updated list of paths to the merged images.
        """
        # Create the main figure and axis
        fig = pylab.figure()
        ax = fig.add_subplot(111)

        # Data for the x-axis and random data for demonstration
        x = range(10)
        data_points = [pylab.randn(10) for _ in range(len(self.split_data))]

        # Plot lines on the same axis and collect them into a list
        lines = []
        filtered_split_data = [
            entry for entry in self.split_data if "FRAME" not in entry
        ]
        for i, data in enumerate(filtered_split_data):
            y = data_points[i]
            label = data.split()[-1]
            type = data.split()[-2]
            if label == "hydrophobic":
                (line,) = ax.plot(
                    x, y, label=data, color=(1.0, 1.0, 0.0), linewidth=5.0
                )  # yellow
            elif label == "hbond":
                if type == "Acceptor":
                    (line,) = ax.plot(
                        x, y, label=data, color=(1.0, 0.6, 0.6), linewidth=5.0
                    )  # light red / pink
                elif type == "Donor":
                    (line,) = ax.plot(
                        x, y, label=data, color=(0.3, 0.5, 1.0), linewidth=5.0
                    )  # light blue
            elif label == "halogen":
                (line,) = ax.plot(
                    x, y, label=data, color=(1.0, 0.0, 0.9), linewidth=5.0
                )  # magenta / hot pink
            elif label == "pistacking":
                (line,) = ax.plot(
                    x, y, label=data, color=(0.0, 0.0, 1.0), linewidth=5.0
                )  # blue
            elif label == "pication":
                (line,) = ax.plot(
                    x, y, label=data, color=(0.0, 0.0, 1.0), linewidth=5.0
                )  # blue
            elif label == "waterbridge":
                (line,) = ax.plot(
                    x, y, label=data, color=(0.0, 1.0, 0.9), linewidth=5.0
                )  # cyan / aqua
            elif label == "metal":
                (line,) = ax.plot(
                    x, y, label=data, color=(1.0, 0.6, 0.0), linewidth=5.0
                )  # orange
            elif label == "saltbridge":
                if type == "NI":
                    (line,) = ax.plot(
                        x, y, label=data, color=(1.0, 0.6, 0.0), linewidth=5.0
                    )  # orange
                elif type == "PI":
                    (line,) = ax.plot(
                        x, y, label=data, color=(0.3, 0.9, 0.8), linewidth=5.0
                    )  # turquoise / teal
            else:
                (line,) = ax.plot(x, y, label=data)
            lines.append(line)

        # Create a separate figure for the legend
        figlegend = pylab.figure(figsize=(8, 6))

        # Add a legend to the subplot (ax) using the lines and full entries as labels
        legend = figlegend.legend(lines, filtered_split_data, loc="center")

        # Set the linewidth of the legend lines to be thicker
        for line in legend.get_lines():
            line.set_linewidth(5.0)

        # Add text above the legend
        figlegend.text(
            0.5, 0.9, f"{self.binding_mode}", ha="center", fontsize=12, weight="bold"
        )
        figlegend.text(
            0.5,
            0.85,
            f"Occurrence {self.occurrence_percent}%",
            ha="center",
            fontsize=12,
            weight="bold",
        )

        # Save the legend figure to a file
        legend_filename = f"{self.binding_mode}_legend.png"
        figlegend.savefig(legend_filename)

        # Read the two images
        image1 = Image.open(f"{self.binding_mode}.png")
        image2 = Image.open(legend_filename)

        # Resize the first image
        image1_size = image1.size
        image2_size = image2.size
        total_width = image1_size[0] + image2_size[0]
        new_image = Image.new("RGB", (total_width, image1_size[1]))
        new_image.paste(image1, (0, 0))
        new_image.paste(image2, (image1_size[0], 0))

        # Save the merged image
        merged_image_filename = f"{self.binding_mode}_merged.png"
        new_image.save(merged_image_filename, "PNG")

        # Append the merged image path to the list
        self.merged_image_paths.append(merged_image_filename)

        # Remove the original files
        os.remove(f"{self.binding_mode}.png")
        os.remove(legend_filename)
        os.remove(f"{self.binding_mode}.svg")

        return self.merged_image_paths


class FigureArranger:
    """
    Arranges multiple merged binding mode figures into a single image.

    Attributes
    ----------
    merged_image_paths : list of str
        List of file paths to pre merged figures.
    output_path : str
        File path where the final arranged figure should be saved.
    """
    def __init__(self, merged_image_paths, output_path):
        self.merged_image_paths = merged_image_paths
        self.output_path = output_path

    def arranged_figure_generation(self):
        """
        Generate an arranged figure by arranging merged images in rows and columns.

        Returns
        -------
        None
        """
        # Open the list of images
        merged_images = [Image.open(path) for path in self.merged_image_paths]

        # Calculate the maximum width and height for the images
        max_width = max(image.size[0] for image in merged_images)
        max_height = max(image.size[1] for image in merged_images)

        # Determine the number of images per row (in your case, 2 images per row)
        images_per_row = 2

        # Calculate the number of rows and columns required
        num_rows = (len(merged_images) + images_per_row - 1) // images_per_row
        total_width = max_width * images_per_row
        total_height = max_height * num_rows

        # Create a new image with the calculated width and height
        big_figure = Image.new(
            "RGB", (total_width, total_height), (255, 255, 255)
        )  # Set background to white

        x_offset = 0
        y_offset = 0

        for image in merged_images:
            # Paste the image onto the big figure
            big_figure.paste(image, (x_offset, y_offset))

            # Update offsets
            x_offset += max_width

            # Move to the next row if necessary
            if x_offset >= total_width:
                x_offset = 0
                y_offset += max_height

        # Save the big figure
        big_figure.save(self.output_path, "PNG")

        # Ensure target directory exists
        target_dir = "Binding_Modes_Markov_States"
        os.makedirs(target_dir, exist_ok=True)

        # Move the file
        new_path = os.path.join(target_dir, os.path.basename(self.output_path))
        os.rename(self.output_path, new_path)

        # Remove the individual image files
        for path in self.merged_image_paths:
            os.remove(path)
