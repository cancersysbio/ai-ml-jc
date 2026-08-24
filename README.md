# AI/ML Journal Club

A journal club exploring AI/ML methods in biology and beyond, with a focus on deep understanding through implementation (if time permits).

## How do we do it 

Understanding an method has several layers:

1. **Just reading the paper** — you may feel like you understand it, but most likely have missed something important.
2. **Implementing the method yourself** — ensures you actually got it and likely provides insights into the implementation tricks that papers gloss over.

<u>Presenters are encouraged to take the second route</u>: implement the method and share a notebook or repository that others can use. Other participants are also encouraged to attempt their own implementations and share the insights during the meeting. 
If you don't have time, feel free to ask to reschedule or at minimum share notes or a presentation with the paper's breakdown. 

## Schedule

- **Frequency:** Every 3 weeks, Mondays at 10:00 AM
- **Start date:** April 20, 2026

| # | Date | Presenter | Paper | Session |
|---|------|-----------|-------|---------|
| 1 | 2026-04-22 | @tulerpetontidae (Artem) | Neural Ordinary Differential Equations (Chen et al., 2018) | [session](sessions/01_neural_odes/) |
| 2 | TBD | TBD | TBD | |
| 3 | TBD | TBD | TBD | |
| 4 | TBD | TBD | TBD | |
| 5 | 2026-07-27 | TBD | TBD | |
| 6 | 2026-08-24 | Maryam Pourmaleki | CytoSignal (Liu et al., *Nature Genetics* 2026) | [session](sessions/06_cytosignal/) |

## Repo Structure

Each session gets its own folder under `sessions/`:

```
sessions/
  01_neural_odes/
    README.md              # paper info, presenter, key concepts
    paper.pdf              # the paper
    implementation.ipynb   # implementation notebook
```

## How to Add a Session

1. Create a new folder under `sessions/` with the format `NN_short_name/`
2. Add the paper PDF as `paper.pdf`
3. Create a `README.md` with paper metadata (title, authors, link, key concepts)
4. Add your implementation notebook as `implementation.ipynb`
5. Update the schedule table in this README

