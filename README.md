# Gnome-Shell Extensions

You will find here some work related to gnome-shell especially the following extensions:

 * French Revolutionary Calendar : a utility that displays the French Republican Calendar.



# French Revolutionary Calendar

This utility displays the current date of the [French Republican Calendar](http://en.wikipedia.org/wiki/French_Republican_Calendar) in the top panel using the [equinox rule](http://en.wikipedia.org/wiki/French_Republican_Calendar#Converting_from_the_Gregorian_Calendar).

It also displays (when clicked) details such as the name of the day and the aspect celebrated this day (often a plant or an instrument of labor).

![Version 8 in action](FRC@jcdubacq.dubacq.fr/screenshot.jpg?raw=true "Version 8 in action")


## Evolution of the extension :

 1. Integration with the real `datetime` box would be really great, but it may involve a lot of code redundancy (fragile) with the original code. The date would appear in the `datetime` button, and the extended info in the date and time panel (below the calendar).
 2. The astronomical computations routines come from [fourmilab](https://www.fourmilab.ch/documents/calendar/). There are plenty of other fun calendars there. The date could be also displayed in Hebrew, Mayan, Islamic or Persion version! So cool... (ongoing ; hebrew and islamic calendars already in).

## Bugs

~~When the date disappears (screen lock), the FRC does not. When the date reappears, the FRC is put to the left instead of to the right.~~ this is no more the case.

The text is not translated, but should it? This is, after all, the **French** Republican Calendar.