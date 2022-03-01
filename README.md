# A personal unpeer-reviewed guide to quantitative FRET

## Brief intro
I wrote this small guide during my Ph.D. while I was developing some Python code to spectrally correct FRET data I measured in order to obtain accurate and objective FRET efficiencies from fluorescence images taking into account the optics used during imaging.

## Relevance
Forster Resonance Energy Transfer (FRET) imaging is often used in biology to dynamically follow biochemical processes within the native context. Lots of new FRET sensors are constantly developed and improved for many different targets. In "traditional" FRET imaging, a 
donor (donor-excitation, donor-emission) and a FRET channel (donor-excitation, acceptor-emission) are acquired and after fluorescence intensity extraction, the FRET "response" is usually expressed as the FRET/Donor ratio. This ratio is mostly not a proxy for "real" FRET efficiency, which actually expressed the physical energy transfer between donor and acceptor. However, this efficiency can be sometimes calculated, either by performing particular experiments with special microscopes, or by "correcting" data obtained with more standardly-available widefield microscopes. I discuss the latter in this small letter.

## Description
The goal of this small letter is to guide the reader (mostly me) through the challenges and warnings when trying to obtaining quantitative FRET efficiencies from your data.

## Software
This letter was written in [TeXMacs](https://www.texmacs.org/tmweb/home/welcome.en.html), mostly because I wanted to explore the tool and use it for something useful :)

## Disclaimer
This is a personal letter, which I mostly wrote to remind myself of the train of thoughts and few calculations I did during my Ph.D. This document has not been peer-reviewed, so it might contain mistakes or incorrect statements.
