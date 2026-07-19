#!/usr/bin/env python3
"""
GIF generator for cellify terminal simulation.
Simulates a cellify session:
1. Querying Materials Project for TiO2
2. Selecting Rutile structure (mp-2657)
3. Converting to conventional cell and building a 2x2x4 supercell.
Outputs to docs/images/demo.gif.
"""

import os
from typing import List, Tuple
from PIL import Image, ImageDraw, ImageFont

# Define path configurations
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
ROOT_DIR = os.path.abspath(os.path.join(SCRIPT_DIR, ".."))
OUTPUT_DIR = os.path.join(ROOT_DIR, "docs/images")
OUTPUT_GIF_PATH = os.path.join(OUTPUT_DIR, "demo.gif")

# Fonts list to try in order
FONT_PATHS = [
    "/System/Library/Fonts/Monaco.ttf",
    "/System/Library/Fonts/Supplemental/Courier New.ttf",
    "/System/Library/Fonts/Supplemental/Andale Mono.ttf",
]

# Color Palette (One Dark Theme)
COLOR_BG_OUTER = "#121212"      # Deep charcoal desktop background
COLOR_BG_INNER = "#1e1e1e"      # Terminal background
COLOR_BG_TITLE = "#2d2d2d"      # Title bar background
COLOR_BORDER = "#333333"        # Window border color
COLOR_TEXT_TITLE = "#abb2bf"    # Title bar text

# macOS Button colors
COLOR_BTN_CLOSE = "#ff5f56"
COLOR_BTN_MIN = "#ffbd2e"
COLOR_BTN_MAX = "#27c93f"

# ANSI Terminal colors
COLOR_PROMPT = "#56b6c2"        # Cyan for prompt
COLOR_DEFAULT = "#abb2bf"       # Light gray for standard output
COLOR_TYPED = "#ffffff"         # White for typed text
COLOR_HIGHLIGHT = "#61afef"     # Blue for important terms
COLOR_SUCCESS = "#98c379"       # Green for success
COLOR_WARNING = "#e5c07b"       # Yellow for warnings

# Dimensions
WIDTH = 800
HEIGHT = 500
MARGIN_X = 10
MARGIN_Y = 10
TITLE_BAR_HEIGHT = 30
PAD_LEFT = 15
PAD_TOP = 15
LINE_HEIGHT = 19
MAX_LINES = 22

# Type definition for terminal segment: (text, color)
Segment = Tuple[str, str]
Line = List[Segment]


def load_font(size: int) -> ImageFont.ImageFont:
    """Loads a monospace font with fallback options.

    If none of the system fonts are found, falls back to default.
    """
    for path in FONT_PATHS:
        if os.path.exists(path):
            try:
                return ImageFont.truetype(path, size)
            except Exception:
                pass
    return ImageFont.load_default()


def draw_terminal(lines: List[Line], cursor_visible: bool) -> Image.Image:
    """Renders the terminal window and its lines to a PIL Image."""
    # 1. Base Image
    img = Image.new("RGB", (WIDTH, HEIGHT), COLOR_BG_OUTER)
    draw = ImageDraw.Draw(img)

    # 2. Outer Window Outline
    draw.rounded_rectangle(
        [(MARGIN_X, MARGIN_Y), (WIDTH - MARGIN_X, HEIGHT - MARGIN_Y)],
        radius=8,
        fill=COLOR_BG_INNER,
        outline=COLOR_BORDER,
        width=1,
    )

    # 3. Title Bar
    draw.rounded_rectangle(
        [(MARGIN_X, MARGIN_Y), (WIDTH - MARGIN_X, MARGIN_Y + TITLE_BAR_HEIGHT)],
        radius=8,
        fill=COLOR_BG_TITLE,
    )
    draw.rectangle(
        [
            (MARGIN_X, MARGIN_Y + TITLE_BAR_HEIGHT // 2),
            (WIDTH - MARGIN_X, MARGIN_Y + TITLE_BAR_HEIGHT),
        ],
        fill=COLOR_BG_TITLE,
    )

    # Re-stroke outline for the top part
    draw.rounded_rectangle(
        [(MARGIN_X, MARGIN_Y), (WIDTH - MARGIN_X, HEIGHT - MARGIN_Y)],
        radius=8,
        outline=COLOR_BORDER,
        width=1,
    )

    # 4. macOS Window Buttons (Close, Minimize, Maximize)
    dot_y = MARGIN_Y + TITLE_BAR_HEIGHT // 2
    dot_radius = 5
    draw.ellipse(
        [
            (25 - dot_radius, dot_y - dot_radius),
            (25 + dot_radius, dot_y + dot_radius),
        ],
        fill=COLOR_BTN_CLOSE,
    )
    draw.ellipse(
        [
            (41 - dot_radius, dot_y - dot_radius),
            (41 + dot_radius, dot_y + dot_radius),
        ],
        fill=COLOR_BTN_MIN,
    )
    draw.ellipse(
        [
            (57 - dot_radius, dot_y - dot_radius),
            (57 + dot_radius, dot_y + dot_radius),
        ],
        fill=COLOR_BTN_MAX,
    )

    # 5. Window Title
    title_text = "cellify — demo"
    font_title = load_font(11)
    title_len = draw.textlength(title_text, font=font_title)
    title_x = WIDTH // 2 - int(title_len) // 2
    title_y = MARGIN_Y + (TITLE_BAR_HEIGHT - 12) // 2
    draw.text(
        (title_x, title_y), title_text, fill=COLOR_TEXT_TITLE, font=font_title
    )

    # 6. Render Lines of Text
    font_terminal = load_font(13)

    visible_lines = lines
    if len(lines) > MAX_LINES:
        visible_lines = lines[-MAX_LINES:]

    start_y = MARGIN_Y + TITLE_BAR_HEIGHT + PAD_TOP
    for i, line in enumerate(visible_lines):
        curr_y = start_y + i * LINE_HEIGHT
        curr_x = MARGIN_X + PAD_LEFT

        for segment_text, color in line:
            draw.text(
                (curr_x, curr_y), segment_text, fill=color, font=font_terminal
            )
            curr_x += int(draw.textlength(segment_text, font=font_terminal))

        # Draw block cursor on the active (last) line
        if cursor_visible and i == len(visible_lines) - 1:
            cursor_height = 14
            cursor_width = 8
            draw.rectangle(
                [
                    (curr_x, curr_y + 1),
                    (curr_x + cursor_width, curr_y + 1 + cursor_height),
                ],
                fill=COLOR_PROMPT,
            )

    return img


def type_text(
    lines: List[Line],
    text: str,
    color: str,
    frames_list: List[Image.Image],
    durations_list: List[int],
    typing_speed: int = 60,
) -> None:
    """Simulates character-by-character typing into the active terminal line."""
    for char in text:
        last_line = lines[-1]
        if not last_line or last_line[-1][1] != color:
            last_line.append((char, color))
        else:
            prev_text, prev_color = last_line[-1]
            last_line[-1] = (prev_text + char, prev_color)

        img = draw_terminal(lines, cursor_visible=True)
        frames_list.append(img)
        durations_list.append(typing_speed)


def blink_cursor(
    lines: List[Line],
    count: int,
    duration: int,
    frames_list: List[Image.Image],
    durations_list: List[int],
) -> None:
    """Generates frames to blink the cursor in a idle state."""
    for _ in range(count):
        # Cursor ON
        img = draw_terminal(lines, cursor_visible=True)
        frames_list.append(img)
        durations_list.append(duration)
        # Cursor OFF
        img = draw_terminal(lines, cursor_visible=False)
        frames_list.append(img)
        durations_list.append(duration)


def print_line(
    lines: List[Line],
    segments: Line,
    delay: int,
    frames_list: List[Image.Image],
    durations_list: List[int],
) -> None:
    """Appends a new line of output to the terminal and records a frame."""
    lines.append(segments)
    img = draw_terminal(lines, cursor_visible=False)
    frames_list.append(img)
    durations_list.append(delay)


def main() -> None:
    """Main sequence builder for generating the demo GIF."""
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    frames: List[Image.Image] = []
    durations: List[int] = []

    # Initial state: Empty prompt
    lines: List[Line] = [[("~ % ", COLOR_PROMPT)]]

    # 1. Initial idle cursor blinking (1.5 seconds)
    blink_cursor(
        lines,
        count=2,
        duration=500,
        frames_list=frames,
        durations_list=durations,
    )

    # 2. Type out the first command
    command_str = (
        "cellify -i TiO2 --conventional -d 2 2 4 -o rutile_supercell.POSCAR"
    )
    type_text(lines, command_str, COLOR_TYPED, frames, durations, typing_speed=65)

    # 3. Wait a moment before running (cursor ON)
    blink_cursor(
        lines,
        count=1,
        duration=400,
        frames_list=frames,
        durations_list=durations,
    )

    # 4. Hitting Enter: network query simulation
    print_line(
        lines,
        [
            (
                "Querying Materials Project OPTIMADE for 'TiO2'...",
                COLOR_DEFAULT,
            )
        ],
        400,
        frames,
        durations,
    )
    print_line(
        lines,
        [("Found 4 structures in Materials Project:", COLOR_DEFAULT)],
        200,
        frames,
        durations,
    )

    # Display structures returned by Materials Project
    mp_structures = [
        [
            ("  [1] ", COLOR_PROMPT),
            ("ID: ", COLOR_DEFAULT),
            ("mp-2657", COLOR_HIGHLIGHT),
            (" | P42/mnm | Vol: 62.4 A^3 | ", COLOR_DEFAULT),
            ("[Stable ★]", COLOR_SUCCESS),
            (" | Rutile", COLOR_TYPED),
        ],
        [
            ("  [2] ", COLOR_PROMPT),
            ("ID: ", COLOR_DEFAULT),
            ("mp-390 ", COLOR_HIGHLIGHT),
            (
                " | I41/amd  | Vol: 136.3 A^3 | E_hull: 0.01 eV | Anatase",
                COLOR_DEFAULT,
            ),
        ],
        [
            ("  [3] ", COLOR_PROMPT),
            ("ID: ", COLOR_DEFAULT),
            ("mp-2658", COLOR_HIGHLIGHT),
            (
                " | A2/m     | Vol: 140.8 A^3 | E_hull: 0.01 eV | Bronze",
                COLOR_DEFAULT,
            ),
        ],
        [
            ("  [4] ", COLOR_PROMPT),
            ("ID: ", COLOR_DEFAULT),
            ("mp-2659", COLOR_HIGHLIGHT),
            (
                " | Pbca     | Vol: 256.1 A^3 | E_hull: 0.02 eV | Brookite",
                COLOR_DEFAULT,
            ),
        ],
    ]
    for struct in mp_structures:
        print_line(lines, struct, 100, frames, durations)

    # Querying COD
    print_line(lines, [("", COLOR_DEFAULT)], 100, frames, durations)
    print_line(
        lines,
        [
            (
                "Querying Crystallography Open Database (COD) OPTIMADE for 'TiO2'...",
                COLOR_DEFAULT,
            )
        ],
        400,
        frames,
        durations,
    )
    print_line(
        lines,
        [
            (
                "Found 40 structures in Crystallography Open Database (COD):",
                COLOR_DEFAULT,
            )
        ],
        200,
        frames,
        durations,
    )

    cod_structures = [
        [
            ("  [5] ", COLOR_PROMPT),
            ("ID: ", COLOR_DEFAULT),
            ("1544349", COLOR_HIGHLIGHT),
            (
                " | P b n m  | Vol: 137.2 A^3 | LixTiO2 x=0.0 ramsdellite-type",
                COLOR_DEFAULT,
            ),
        ],
        [
            ("  [6] ", COLOR_PROMPT),
            ("ID: ", COLOR_DEFAULT),
            ("1557789", COLOR_HIGHLIGHT),
            (
                " | P 1 2/c 1| Vol: 121.5 A^3 | Riesite IMA-2015-110",
                COLOR_DEFAULT,
            ),
        ],
    ]
    for struct in cod_structures:
        print_line(lines, struct, 100, frames, durations)

    print_line(
        lines,
        [
            (
                "  Warning: Only the top 5 most relevant structures are shown. "
                "There are 35 more structures in Crystallography Open Database (COD).",
                COLOR_WARNING,
            )
        ],
        200,
        frames,
        durations,
    )
    print_line(lines, [("", COLOR_DEFAULT)], 150, frames, durations)

    # 5. Show select prompt
    select_prompt_line = [("Select a structure (1-10): ", COLOR_DEFAULT)]
    print_line(lines, select_prompt_line, 400, frames, durations)

    # Idle cursor on prompt
    blink_cursor(
        lines,
        count=1,
        duration=400,
        frames_list=frames,
        durations_list=durations,
    )

    # Type choice '1'
    type_text(lines, "1", COLOR_TYPED, frames, durations, typing_speed=150)

    # Pause before execution
    blink_cursor(
        lines,
        count=1,
        duration=300,
        frames_list=frames,
        durations_list=durations,
    )

    # 6. Press Enter: Download and process Rutile
    print_line(
        lines,
        [
            (
                "Downloading structure from Materials Project (ID: mp-2657)...",
                COLOR_DEFAULT,
            )
        ],
        500,
        frames,
        durations,
    )
    print_line(
        lines,
        [("  Formula: ", COLOR_DEFAULT), ("TiO2", COLOR_TYPED)],
        150,
        frames,
        durations,
    )
    print_line(
        lines,
        [("  Volume:  ", COLOR_DEFAULT), ("62.433 A^3", COLOR_HIGHLIGHT)],
        150,
        frames,
        durations,
    )
    print_line(
        lines,
        [("  Number of atoms: ", COLOR_DEFAULT), ("6", COLOR_HIGHLIGHT)],
        250,
        frames,
        durations,
    )
    print_line(
        lines,
        [
            (
                "Converting structure to standard conventional cell...",
                COLOR_DEFAULT,
            )
        ],
        400,
        frames,
        durations,
    )
    print_line(
        lines,
        [
            (
                "Generating supercell with diagonal scaling: ",
                COLOR_DEFAULT,
            ),
            ("[2, 2, 4]", COLOR_HIGHLIGHT),
        ],
        500,
        frames,
        durations,
    )

    # Final summary
    print_line(
        lines, [("Final structure summary:", COLOR_SUCCESS)], 150, frames, durations
    )
    print_line(
        lines,
        [("  Formula: ", COLOR_DEFAULT), ("Ti32O64", COLOR_TYPED)],
        100,
        frames,
        durations,
    )
    print_line(
        lines,
        [("  Volume:  ", COLOR_DEFAULT), ("998.928 A^3", COLOR_HIGHLIGHT)],
        100,
        frames,
        durations,
    )
    print_line(
        lines,
        [("  Number of atoms: ", COLOR_DEFAULT), ("96", COLOR_HIGHLIGHT)],
        100,
        frames,
        durations,
    )
    print_line(
        lines, [("  Lattice constants:", COLOR_DEFAULT)], 100, frames, durations
    )
    print_line(
        lines,
        [
            (
                "    a = 9.1874 A, b = 9.1874 A, c = 11.8348 A",
                COLOR_DEFAULT,
            )
        ],
        100,
        frames,
        durations,
    )
    print_line(
        lines,
        [
            (
                "    alpha = 90.00 deg, beta = 90.00 deg, gamma = 90.00 deg",
                COLOR_DEFAULT,
            )
        ],
        150,
        frames,
        durations,
    )
    print_line(lines, [("", COLOR_DEFAULT)], 100, frames, durations)
    print_line(
        lines,
        [
            ("Saving final structure to: ", COLOR_DEFAULT),
            ("rutile_supercell.POSCAR", COLOR_HIGHLIGHT),
        ],
        300,
        frames,
        durations,
    )
    print_line(lines, [("Success!", COLOR_SUCCESS)], 1000, frames, durations)

    # Next prompt
    print_line(lines, [("~ % ", COLOR_PROMPT)], 200, frames, durations)

    # Final blinks on the terminal before restarting loop (4 seconds pause)
    blink_cursor(
        lines,
        count=4,
        duration=500,
        frames_list=frames,
        durations_list=durations,
    )

    # 7. Compile and save the GIF
    print(f"Compiling {len(frames)} frames to: {OUTPUT_GIF_PATH}")
    frames[0].save(
        OUTPUT_GIF_PATH,
        save_all=True,
        append_images=frames[1:],
        duration=durations,
        loop=0,
        optimize=True,
    )
    print("Demo GIF successfully saved!")


if __name__ == "__main__":
    main()
