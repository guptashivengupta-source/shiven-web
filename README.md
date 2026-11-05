# Shiven Gupta website

Personal website for Shiven Gupta, senior at Prosper High School (Class of 2027) — UNT IDEA Lab student researcher, two-time AIME qualifier, All-State tenor and club basketball player.

Research code: https://github.com/guptashivengupta-source/genie3_brca_project

## Files

- index.html - main page
- resume.html - resume page
- css/style.css - all the styles
- js/main.js - menu, punnett square, slider, prime spiral, copy button
- python/make_figures.py - makes the research charts
- data/ - chart numbers and the free throw log
- images/ - photo and the tab icon

## How to open

Double click index.html. Or run this in the folder and go to http://localhost:8000

```
python -m http.server 8000
```

## Update the charts

```
python python/make_figures.py
```

## Put it online

Drag the folder into https://app.netlify.com/drop, or upload it to GitHub Pages.
