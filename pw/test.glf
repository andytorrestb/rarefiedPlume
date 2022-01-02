#############################################################################
#
# (C) 2021 Cadence Design Systems, Inc. All rights reserved worldwide.
#
# This sample script is not supported by Cadence Design Systems, Inc.
# It is provided freely for demonstration purposes only.
# SEE THE WARRANTY DISCLAIMER AT THE BOTTOM OF THIS FILE.
#
#############################################################################

# --------------------------------------------------------------------------
#
# ConCircleCenterAndRaduis.glf
#
# Script with Tk interface to create connectors for a circle given the
# center point and radius. The circle consists of two connectors, each
# of 180 degrees of arc. The circle is created in the current drawing
# plane or in one of the three Cartesian planes. An optional dimension
# can be specified.
#
# --------------------------------------------------------------------------

package require PWI_Glyph 2

# enable Tk in Glyph 2
pw::Script loadTk

set opt(Center)    [pwu::Vector3 set 0 0 0]
set opt(Radius)     1
set opt(Dimension) ""
set opt(Plane)     view

if { ! [pw::Application isInteractive] } {
  set opt(Plane) xy
}

set color(Valid) SystemWindow
set color(Invalid) MistyRose

# --------------------------------------------------------------------------
# Validate and convert a list of 3 values to a Vector3

proc toVec3 { coord } {
  set flot {([+-]?\d*\.?\d+(?:e[+-]?\d+)?)}
  set sep {(?:(?:\s*,\s*)|(?:\s+))}
  set re "^\\s*${flot}${sep}${flot}${sep}${flot}\\s*$"
  if { [catch { regexp $re $coord a x y z } rc] } {
    return -code error "regexp error: $rc"
  } elseif { ! $rc } {
    return -code error "no match"
  } elseif { [catch { expr \
      [string is double -strict $x] && \
      [string is double -strict $y] && \
      [string is double -strict $z] } rc] } {
    return -code error "not reals"
  }
  return [pwu::Vector3 set $x $y $z]
}

# --------------------------------------------------------------------------
# Validate a coordinate input field

proc validateVec3 { w text action } {
  if { $action == -1 } { return 1 }

  if { ! [catch { toVec3 $text }] } {
    $w configure -bg $::color(Valid)
  } else {
    $w configure -bg $::color(Invalid)
  }

  updateButtons
  return 1
}

# --------------------------------------------------------------------------
# Validate radius input field

proc validateRadius { w text action } {
  if { $action != -1 } {
    if { ! [catch { expr double($text) } text] && 0.0 < $text } {
      $w configure -bg $::color(Valid)
    } else {
      $w configure -bg $::color(Invalid)
    }
    updateButtons
  }
  return 1
}

# --------------------------------------------------------------------------
# Validate dimension input field

proc validateDim { w text action } {
  if { $action != -1 } {
    if { [string length $text] == 0 || ( ! [catch { expr int($text) } text] && \
        (0 == $text || 4 <= $text)) } {
      $w configure -bg $::color(Valid)
    } else {
      $w configure -bg $::color(Invalid)
    }
    updateButtons
  }
  return 1
}

#---------------------------------------------------------------------------
# Create the connectors and dimension them

proc makeConnector { } {
  global opt

  if { [catch { toVec3 $opt(Center) } c] } {
    focus .inputs.entcoord
    return 0
  }
  if { [catch { expr double($opt(Radius)) } r] || $r <= 0.0 } {
    focus .inputs.entrad
    return 0
  }
  set dim $opt(Dimension)
  if { [string length $opt(Dimension)] == 0 } { set dim 0 }
  if { [catch { expr int($dim) } dim] || ($dim < 4 && $dim != 0) } {
    focus .inputs.entdim
    return 0
  }

  set point(0) [pwu::Vector3 set $r 0.0 0.0]
  set point(1) [pwu::Vector3 set  [expr 0.717* $r] [expr 0.717* $r] 0.0]
  set point(2) [pwu::Vector3 set 0.0 $r 0.0]

  # create the semi-circle connectors in the the selected plane
  switch $opt(Plane) {
    view {
      # current plane of view
      set view [pw::Display getCurrentView]
      set axis [lindex $view 2]
      set angle [lindex $view 3]
    }
    xy {
      # XY plane
      set axis { 0 0 1 }
      set angle 0.0
    }
    yz {
      # YZ plane
      set axis { 0 1 0 }
      set angle 90.0
    }
    xz {
      # XZ plane
      set axis { 1 0 0 }
      set angle 90.0
    }
  }

  set xform [pwu::Transform inverse [pwu::Transform rotation $axis $angle]]

  set point(0) [pwu::Vector3 add [pwu::Transform apply $xform $point(0)] $c]
  set point(1) [pwu::Vector3 add [pwu::Transform apply $xform $point(1)] $c]
  set point(2) [pwu::Vector3 add [pwu::Transform apply $xform $point(2)] $c]
  
  set mode [pw::Application begin Create]

  if { [catch {
    set con1 [pw::Connector create]
    set seg [pw::SegmentCircle create]
    $seg addPoint $point(0)
    $seg addPoint $point(2)
    $seg setShoulderPoint $point(1)
    $con1 addSegment $seg

    if { $dim >= 4 } {
      set dim1 [expr {1 + ($dim)/2}]
      set dim2 [expr {$dim + 2 - $dim1}]
      $con1 setDimension $dim1
    }
  } msg] } {
    tk_messageBox -icon error -title "Could not create connectors" \
      -message $msg -type ok
    $mode abort
    return 0
  }
  $mode end
  return 1
}

# --------------------------------------------------------------------------
# Enable/disable action buttons based on input validity

proc updateButtons { } {
  set ok [expr [string equal [.inputs.entcoord cget -bg] $::color(Valid)] && \
    [string equal [.inputs.entrad cget -bg] $::color(Valid)] && \
    [string equal [.inputs.entdim cget -bg] $::color(Valid)]]
    
  if { $ok } {
    .buttons.ok configure -state normal
  } else {
    .buttons.ok configure -state disabled
  }
}

# --------------------------------------------------------------------------
# Generate and lay out field/label widget pair

proc makeInputField { parent name title varname {width 7} {valid ""}} {
  set l $parent.lbl$name
  set e $parent.ent$name
  ttk::label $l -text $title
  entry $e -textvariable $varname -width $width
  if { ! [string equal $valid ""] } {
    $e configure -validate all
    $e configure -validatecommand $valid
  }
  grid $l $e
  grid configure $l -sticky e -padx 2
  grid configure $e -sticky w -padx 4

  return $parent.$name
}

# --------------------------------------------------------------------------
# Build the interface

proc makeWindow { } {
  global opt

  wm title . "Circular Connector"

  label .title -text "Circular Connector"
  set font [.title cget -font]
  .title configure -font \
    [font create -family [font actual $font -family] -weight bold]
  pack .title -expand 1 -side top

  pack [frame .hr1 -relief sunken -height 2 -bd 1] -side top -padx 2 \
    -fill x -pady 1
  pack [frame .inputs] -fill x -padx 2 -expand 1

  makeInputField .inputs coord "Center Coordinate:" opt(Center) 20 \
    "validateVec3 %W %P %d"
  makeInputField .inputs rad "Radius:" opt(Radius) 20 \
    "validateRadius %W %P %d"
  makeInputField .inputs dim "Dimension:" opt(Dimension) 20 \
    "validateDim %W %P %d"

  ttk::labelframe .inputs.p -text "Creation Plane"
  ttk::radiobutton .inputs.p.v -text "View plane" -variable opt(Plane) \
    -value view
  ttk::radiobutton .inputs.p.xy -text "XY" -variable opt(Plane) -value xy
  ttk::radiobutton .inputs.p.yz -text "YZ" -variable opt(Plane) -value yz
  ttk::radiobutton .inputs.p.xz -text "XZ" -variable opt(Plane) -value xz

  if { ! [pw::Application isInteractive] } {
    .inputs.p.v configure -state disabled
  }

  grid .inputs.p -columnspan 2 -sticky ew -pady 6 -padx 12
  grid .inputs.p.v .inputs.p.xy -sticky w
  grid .inputs.p.yz .inputs.p.xz -sticky w

  grid columnconfigure .inputs 1 -weight 1

  pack [frame .hr2 -relief sunken -height 2 -bd 1] \
    -side top -padx 2 -fill x -pady 1

  pack [frame .buttons] -fill x -padx 2 -pady 1
  pack [button .buttons.cancel -text "Cancel" -command { exit }] \
    -side right -padx 2
  pack [button .buttons.ok -text "OK" -command { if [makeConnector] exit }] \
    -side right -padx 2
  pack [label .buttons.logo -image [cadenceLogo] -bd 0 -relief flat] \
    -side left -padx 5

  updateButtons

  bind . <Control-Return> { .buttons.ok invoke }
  bind . <Escape> { .buttons.cancel invoke }
}

# --------------------------------------------------------------------------
# Cadence Design Systems logo image

proc cadenceLogo {} {
  set logoData "
R0lGODlhgAAYAPQfAI6MjDEtLlFOT8jHx7e2tv39/RYSE/Pz8+Tj46qoqHl3d+vq62ZjY/n4+NT
T0+gXJ/BhbN3d3fzk5vrJzR4aG3Fubz88PVxZWp2cnIOBgiIeH769vtjX2MLBwSMfIP///yH5BA
EAAB8AIf8LeG1wIGRhdGF4bXD/P3hwYWNrZXQgYmVnaW49Iu+7vyIgaWQ9Ilc1TTBNcENlaGlIe
nJlU3pOVGN6a2M5ZCI/PiA8eDp4bXBtdGEgeG1sbnM6eD0iYWRvYmU6bnM6bWV0YS8iIHg6eG1w
dGs9IkFkb2JlIFhNUCBDb3JlIDUuMC1jMDYxIDY0LjE0MDk0OSwgMjAxMC8xMi8wNy0xMDo1Nzo
wMSAgICAgICAgIj48cmRmOlJERiB4bWxuczpyZGY9Imh0dHA6Ly93d3cudy5vcmcvMTk5OS8wMi
8yMi1yZGYtc3ludGF4LW5zIyI+IDxyZGY6RGVzY3JpcHRpb24gcmY6YWJvdXQ9IiIg/3htbG5zO
nhtcE1NPSJodHRwOi8vbnMuYWRvYmUuY29tL3hhcC8xLjAvbW0vIiB4bWxuczpzdFJlZj0iaHR0
cDovL25zLmFkb2JlLmNvbS94YXAvMS4wL3NUcGUvUmVzb3VyY2VSZWYjIiB4bWxuczp4bXA9Imh
0dHA6Ly9ucy5hZG9iZS5jb20veGFwLzEuMC8iIHhtcE1NOk9yaWdpbmFsRG9jdW1lbnRJRD0idX
VpZDoxMEJEMkEwOThFODExMUREQTBBQzhBN0JCMEIxNUM4NyB4bXBNTTpEb2N1bWVudElEPSJ4b
XAuZGlkOkIxQjg3MzdFOEI4MTFFQjhEMv81ODVDQTZCRURDQzZBIiB4bXBNTTpJbnN0YW5jZUlE
PSJ4bXAuaWQ6QjFCODczNkZFOEI4MTFFQjhEMjU4NUNBNkJFRENDNkEiIHhtcDpDcmVhdG9yVG9
vbD0iQWRvYmUgSWxsdXN0cmF0b3IgQ0MgMjMuMSAoTWFjaW50b3NoKSI+IDx4bXBNTTpEZXJpZW
RGcm9tIHN0UmVmOmluc3RhbmNlSUQ9InhtcC5paWQ6MGE1NjBhMzgtOTJiMi00MjdmLWE4ZmQtM
jQ0NjMzNmNjMWI0IiBzdFJlZjpkb2N1bWVudElEPSJ4bXAuZGlkOjBhNTYwYTM4LTkyYjItNDL/
N2YtYThkLTI0NDYzMzZjYzFiNCIvPiA8L3JkZjpEZXNjcmlwdGlvbj4gPC9yZGY6UkRGPiA8L3g
6eG1wbWV0YT4gPD94cGFja2V0IGVuZD0iciI/PgH//v38+/r5+Pf29fTz8vHw7+7t7Ovp6Ofm5e
Tj4uHg397d3Nva2djX1tXU09LR0M/OzczLysnIx8bFxMPCwcC/vr28u7q5uLe2tbSzsrGwr66tr
KuqqainpqWko6KhoJ+enZybmpmYl5aVlJOSkZCPjo2Mi4qJiIeGhYSDgoGAf359fHt6eXh3dnV0
c3JxcG9ubWxramloZ2ZlZGNiYWBfXl1cW1pZWFdWVlVUU1JRUE9OTUxLSklIR0ZFRENCQUA/Pj0
8Ozo5ODc2NTQzMjEwLy4tLCsqKSgnJiUkIyIhIB8eHRwbGhkYFxYVFBMSERAPDg0MCwoJCAcGBQ
QDAgEAACwAAAAAgAAYAAAF/uAnjmQpTk+qqpLpvnAsz3RdFgOQHPa5/q1a4UAs9I7IZCmCISQwx
wlkSqUGaRsDxbBQer+zhKPSIYCVWQ33zG4PMINc+5j1rOf4ZCHRwSDyNXV3gIQ0BYcmBQ0NRjBD
CwuMhgcIPB0Gdl0xigcNMoegoT2KkpsNB40yDQkWGhoUES57Fga1FAyajhm1Bk2Ygy4RF1seCjw
vAwYBy8wBxjOzHq8OMA4CWwEAqS4LAVoUWwMul7wUah7HsheYrxQBHpkwWeAGagGeLg717eDE6S
4HaPUzYMYFBi211FzYRuJAAAp2AggwIM5ElgwJElyzowAGAUwQL7iCB4wEgnoU/hRgIJnhxUlpA
SxY8ADRQMsXDSxAdHetYIlkNDMAqJngxS47GESZ6DSiwDUNHvDd0KkhQJcIEOMlGkbhJlAK/0a8
NLDhUDdX914A+AWAkaJEOg0U/ZCgXgCGHxbAS4lXxketJcbO/aCgZi4SC34dK9CKoouxFT8cBNz
Q3K2+I/RVxXfAnIE/JTDUBC1k1S/SJATl+ltSxEcKAlJV2ALFBOTMp8f9ihVjLYUKTa8Z6GBCAF
rMN8Y8zPrZYL2oIy5RHrHr1qlOsw0AePwrsj47HFysrYpcBFcF1w8Mk2ti7wUaDRgg1EISNXVwF
lKpdsEAIj9zNAFnW3e4gecCV7Ft/qKTNP0A2Et7AUIj3ysARLDBaC7MRkF+I+x3wzA08SLiTYER
KMJ3BoR3wzUUvLdJAFBtIWIttZEQIwMzfEXNB2PZJ0J1HIrgIQkFILjBkUgSwFuJdnj3i4pEIlg
eY+Bc0AGSRxLg4zsblkcYODiK0KNzUEk1JAkaCkjDbSc+maE5d20i3HY0zDbdh1vQyWNuJkjXnJ
C/HDbCQeTVwOYHKEJJwmR/wlBYi16KMMBOHTnClZpjmpAYUh0GGoyJMxya6KcBlieIj7IsqB0ji
5iwyyu8ZboigKCd2RRVAUTQyBAugToqXDVhwKpUIxzgyoaacILMc5jQEtkIHLCjwQUMkxhnx5I/
seMBta3cKSk7BghQAQMeqMmkY20amA+zHtDiEwl10dRiBcPoacJr0qjx7Ai+yTjQvk31aws92JZ
Q1070mGsSQsS1uYWiJeDrCkGy+CZvnjFEUME7VaFaQAcXCCDyyBYA3NQGIY8ssgU7vqAxjB4EwA
DEIyxggQAsjxDBzRagKtbGaBXclAMMvNNuBaiGAAA7"

  return [image create photo -format GIF -data $logoData]
}

# Create the Tk window and place it
makeWindow
::tk::PlaceWindow . widget

tkwait window .

#############################################################################
#
# This file is licensed under the Cadence Public License Version 1.0 (the
# "License"), a copy of which is found in the included file named "LICENSE",
# and is distributed "AS IS." TO THE MAXIMUM EXTENT PERMITTED BY APPLICABLE
# LAW, CADENCE DISCLAIMS ALL WARRANTIES AND IN NO EVENT SHALL BE LIABLE TO
# ANY PARTY FOR ANY DAMAGES ARISING OUT OF OR RELATING TO USE OF THIS FILE.
# Please see the License for the full text of applicable terms.
#
#############################################################################