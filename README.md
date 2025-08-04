# Reclone Open DNA Collections

> A global collaboration for equitable access to biotechnology

The [Reclone Reagent Collaboration Network](https://reclone.org) has worked for the creation and disponibilization of open DNA collections. While these collections were originally distributed by [Freegenes](https://stanford.freegenes.org)
), now some parts of this project are separating into their own organization structure. The parts in this repository will receive a new ID, and this is where we will centralize the issues board for any changes proposed for the DNA or documentation on these collections.

Here you can find the genbank files and metadata for the plasmids in different collections:
- Open Enzymes Collection
- Open Reporters Collection
- Open Yeast Collection
- Open Plasmids
- E.coli Protein Expression Toolkit
- Molecular Diagnostics Toolkit

## Interactive Database

🧬 **[Launch the Interactive Database](https://your-app-name.streamlit.app)** (deployed on Streamlit Community Cloud)

The interactive database provides:
- **Search & Browse**: Find DNA parts across all collections
- **Analytics**: Visualize collection statistics and sequence properties
- **Part Details**: View detailed information including sequences
- **Data Export**: Download data in various formats

### Local Development

To run the database locally:

```bash
pip install -r requirements.txt
streamlit run streamlit_app.py
```

## Sheets Backend

The backend database for Freegenes lives here: https://docs.google.com/spreadsheets/d/1LZCXzBtgey9xv5OH7YGYgp8UMJ27Eyj1aF9IhAW6M6o/edit?pli=1#gid=954552604

Although it's not currently our intention to use it directly for Reclone, this has important information on every DNA part we have
