"""Place the existing themed X and Z fit panels side by side in one vector PDF."""
from pathlib import Path
from pypdf import PdfReader, PdfWriter, Transformation

STUDY = Path(__file__).resolve().parents[1]
SOURCES = [
    STUDY / 'figures/x-one-step-local-fit.pdf',
    STUDY.parent / 'toric3d-h6-z-one-step-2026-10-15/figures/z-one-step-local-fit.pdf',
]
OUTPUT = STUDY / 'figures/one-step-feedback-x-z.pdf'


def main():
    readers = [PdfReader(path) for path in SOURCES]
    assert all(len(reader.pages) == 1 for reader in readers)
    panels = [reader.pages[0] for reader in readers]
    height = max(float(page.mediabox.height) for page in panels)
    scales = [height / float(page.mediabox.height) for page in panels]
    widths = [float(page.mediabox.width) * scale for page, scale in zip(panels, scales)]
    margin, gap = 8.0, 12.0
    writer = PdfWriter()
    page = writer.add_blank_page(width=sum(widths) + 2 * margin + gap, height=height + 2 * margin)
    x = margin
    for panel, scale, width in zip(panels, scales, widths):
        transform = Transformation().scale(scale).translate(x, margin)
        page.merge_transformed_page(panel, transform, expand=False)
        x += width + gap
    writer.add_metadata({
        '/Title': 'H+6 one-step feedback: X and Z fit panels',
        '/Subject': 'Existing themed panels combined without resampling, refitting, or rasterization.',
    })
    with OUTPUT.open('wb') as stream:
        writer.write(stream)
    result = PdfReader(OUTPUT)
    assert len(result.pages) == 1
    text = result.pages[0].extract_text()
    compact_text = ''.join(text.split())
    for panel in panels:
        for line in panel.extract_text().splitlines():
            # PDF text extraction may insert spaces after page scaling.
            assert ''.join(line.split()) in compact_text, f'Missing source-panel text: {line!r}'
    print('Combined both complete vector panels into one page; all source-panel text retained.')


if __name__ == '__main__':
    main()
