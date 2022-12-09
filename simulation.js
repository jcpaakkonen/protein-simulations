/*
	PROTEIN THERMODYNAMICS SIMULATIONS
	Copyright (C) 2021–2022 Johan Pääkkönen, Juha Rouvinen, University of Eastern Finland
	
	This program is free software; you can redistribute it and/or
	modify it under the terms of the GNU General Public License
	as published by the Free Software Foundation; either version 2
	of the License, or (at your option) any later version.
	
	This program is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	GNU General Public License for more details.
	
	You should have received a copy of the GNU General Public License
	along with this program; if not, write to the Free Software
	Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
*/

/*
	The GNU General Public License, version 2, is available at:
	(1) the file LICENSE in this folder
	(2) https://www.gnu.org/licenses/old-licenses/gpl-2.0.html
*/

// If you use this in your own research, please cite our paper: https://doi.org/10.1021/acsomega.2c00560

"use strict";

// This controls the scales of the slider values.
// Adjust the selector .valuecol (lines 40–49) of style.css accordingly.
// false: logarithmic scale, used in development versions and my PhD thesis
// true: more humanly convenient scale, used in the final version
var newdecadescale = true;

// Set to true to display indicators of the pie segments.
// I made it only for preparing figures for my PhD thesis,
// so do not complain if it behaves weirdly.
var extpiemode = false;

var S_0 = 1;
var E_0 = 1;
var Q_0 = 1;
var K_D = 1;
var K_D2 = 1;

var scale_absolute = 1;
var xscale_alternative = false;

var fwidth_old = 0;
var fheight_old = 0;

var xaxistype = 0;
var yaxistype = 0;
var axistype_lin = 0;
var axistype_log = 1;
var axistype_time = 2;

var curves = [];
var labels = [];
var legends = [];
var colours = [];
var pie = [];
var datapoints = [];

var xmin = 0, xmax = 0;
var ymin = 0, ymax = 0;

var appmode = -1;
var appmode_listonly  = -1;
var appmode_ligands   = 0;
var appmode_receptors = 1;
var appmode_ligand    = 2;
var appmode_homodimer = 3;

var value_digits = newdecadescale ? 1 : 2;

var svg_xmlns = "http://www.w3.org/2000/svg";
var xml_xmlns = "http://www.w3.org/XML/1998/namespace";

var decadeshift = 0;
var extmode = false;

function preinit(mode)
{
	appmode = mode;
	
	var linkarea = document.getElementById("linkarea");
	var credits = document.getElementById("credits");
	var brnode, anode, textnode;
	var first = true;
	var i;
	
	var qm = window.location.href.lastIndexOf("?");
	
	if(qm >= 0)
	{
		var old = qm;
		var amp;
		
		do
		{
			amp = window.location.href.indexOf("&", old + 1);
			var ss = window.location.href.substring(old + 1, amp >= 0 ? amp : undefined);
			if(ss === "ext") extmode = true;
			old = amp;
		}
		while(amp >= 0);
	}
	
	if(extmode)
	{
		var links = document.getElementById("linkarea").getElementsByTagName("a");
		for(i = 0; i < links.length; i++) links[i].href += "?ext";
		
		document.getElementById("ext_extmessage").style.display = "none";
	}
	
	if(!(appmode >= 0)) return;
	
	if(extmode)
	{
		document.getElementById("ext_dsrow").style.display = "table-row";
		document.getElementById("ext_lrrow").style.display = "initial";
		document.getElementById("ext_lrbr").style.display = "initial";
		document.getElementById("ext_pisrow").style.display = "initial";
		document.getElementById("ext_pisbr").style.display = "initial";
	}
	
	var figure_div = document.getElementById("figure_div");
	
	figure_div.style.width = "625px";
	figure_div.style.minWidth = "450px";
	figure_div.style.height = "450px";
	figure_div.style.minHeight = "300px";
	
	/*
	figure_div.style.width = "585px";
	figure_div.style.minWidth = "585px";
	figure_div.style.height = "300px";
	figure_div.style.minHeight = "300px";
	*/
	
	var toggle_calcdiv_checkbox = document.getElementById("toggle_calcdiv_checkbox")
	if(toggle_calcdiv_checkbox)
	{
		// crude hack, but it works...
		// because the function will reverse this back immediately
		toggle_calcdiv_checkbox.checked = !toggle_calcdiv_checkbox.checked;
		toggle_calcdiv();
	}
}

function init()
{
	var i;
	
	if(!(appmode >= 0)) return;
	
	if(extmode) decadeshift_changed(true);
	
	var sliders = document.getElementsByClassName("slider");
	
	for(i = 0; i < sliders.length; i++)
	{
		slider_input(Number(sliders[i].id.substr(6)), true);
	}
	
	var radio1 = document.getElementById("radio1");
	scale_absolute = (radio1 && radio1.checked) ? 1 : 0;
	
	var radio5 = document.getElementById("radio5");
	if(radio5)
	{
		var radio6 = document.getElementById("radio6");
		var radio7 = document.getElementById("radio7");
		
		if(radio5.checked) scale_absolute = 1;
		if(radio6.checked) scale_absolute = 0;
		if(radio7.checked) scale_absolute = 2;
	}
	
	var radio3 = document.getElementById("radio3");
	xscale_alternative = !(radio3 && radio3.checked);
	
	radio_input(xscale_alternative ? 4 : 3, true);
	
	var figure_div = document.getElementById("figure_div");
	fwidth_old = extract_px(figure_div.style.width);
	fheight_old = extract_px(figure_div.style.height);

	document.body.addEventListener("mousemove", mouse_moved);
	document.body.addEventListener("mouseup", mouse_released);
	
	var databox = document.getElementById("databox");
	if(databox) databox.addEventListener("blur", data_changed);
	
	if(extmode)
	{
		document.getElementById("ext_dsinput").addEventListener("input", function(){decadeshift_changed()});
		document.getElementById("ext_dsinput").addEventListener("blur", decadeshift_blur);
	}
	
	var calcbutton = document.getElementById("calcbutton");
	if(calcbutton) calcbutton.addEventListener("click",
		function()
		{
			var label = document.getElementById("calcstatus");
			if(label)
			{
				label.innerHTML = "Calculating, please wait...";
				label.style.color = "blue";
			}
			
			// This is so that the label has time to change
			// before the code enters a heavy calculation.
			window.setTimeout(calculate_lsf, 100);
		}
	);
	
	document.getElementById("export_button").addEventListener("click", export_svg);
	
	data_changed(false);
	update();
}

function valuespan_click(num)
{
	var value_ele = document.getElementById("value" + num);
	var input_ele = document.getElementById("value_input" + num);
	var t;
	
	switch(num)
	{
		case 3: t = S_0; break;
		case 5: t = E_0; break;
		case 7: t = K_D; break;
		case 9: t = K_D2; break;
		case 10: t = Q_0; break;
		default: break;
	}
	
	input_ele.value = t.toExponential(value_digits);
	
	if(t >= 1e-9 && t <= 1e-1)
		input_ele.style.color = "black";
	else if(t >= 1e-18 && t <= 1e3)
		input_ele.style.color = "darkgoldenrod";
	else
		input_ele.style.color = "red";
	
	value_ele.style.display = "none";
	input_ele.style.display = "inline";
	input_ele.focus();
}

function valueinput_update(num)
{
	var input_ele = document.getElementById("value_input" + num);
	var t = Number(input_ele.value);
	
	if(t >= 1e-9 && t <= 1e-1)
		input_ele.style.color = "black";
	else if(t >= 1e-18 && t <= 1e3)
		input_ele.style.color = "darkgoldenrod";
	else
		input_ele.style.color = "red";
	
	if(t >= 1e-18 && t <= 1e3) switch(num)
	{
		case 3: S_0 = t; break;
		case 5: E_0 = t; break;
		case 7: K_D = t; break;
		case 9: K_D2 = t; break;
		case 10: Q_0 = t; break;
		default: break;
	}
	
	slider_input(num, false, true);
}

function valueinput_blur(num, keycode)
{
	var input_ele = document.getElementById("value_input" + num);
	if(input_ele.style.display === "none") return;
	
	if(keycode === undefined || keycode === 13 /* RETURN/ENTER */ || keycode === 27 /* ESC */)
	{
		var value_ele = document.getElementById("value" + num);
		
		//valueinput_update(num);
		
		value_ele.style.display = "inline";
		input_ele.style.display = "none";
		input_ele.blur();
	}
}

function decadeshift_changed(noupdate)
{
	var ext_dsinput = document.getElementById("ext_dsinput");
	var d = Number(ext_dsinput.value);
	var i;
	
	if(ext_dsinput.value !== "" && Number.isInteger(d) && d >= Number(ext_dsinput.min) && d <= Number(ext_dsinput.max))
	{
		ext_dsinput.style.color = "initial";
		decadeshift = d;
		
		if(!noupdate)
		{
			var sliders = document.getElementsByClassName("slider");
			
			for(i = 0; i < sliders.length; i++)
			{
				slider_input(Number(sliders[i].id.substr(6)), true);
			}
			
			update();
		}
	}
	else
	{
		ext_dsinput.style.color = "red";
	}
}

function decadeshift_blur()
{
	var ext_dsinput = document.getElementById("ext_dsinput");
	ext_dsinput.style.color = "initial";
	ext_dsinput.value = decadeshift;
}

function render_text_svg(ele, str)
{
	var fontfamily = "";
	var normalfont = "12px";
	var scriptfont = "10px";
	
	var supbl = -4;
	var subpl = 2.5;
	
	var nextchar = 0;
	var prevchar = -1;
	var search = 0;
	var token = "";
	var baseline = 0;
	var c;
	
	delete_all_children(ele);
	
	do
	{
		search = str.substr(prevchar + 1).search(/[\uEEE0-\uEEEF]/);
		nextchar = prevchar + 1 + search;
		token = str.substring(prevchar + 1, (search === -1) ? undefined : nextchar);
		
		c = (prevchar === -1 ? 0xEEE0 : str.charCodeAt(prevchar));
		
		function new_textspan(bold, italic, script, align, content)
		{
			if(!content) return;

			var span = document.createElementNS(svg_xmlns, "tspan");
			if(bold) span.style.fontWeight = "bold";
			if(italic) span.style.fontStyle = "italic";
			if(script) span.style.fontSize = "10px";
			if(align - baseline !== 0) span.setAttribute("dy", align - baseline);
			baseline = align;

			var textnode = document.createTextNode(content);
			span.appendChild(textnode);

			ele.appendChild(span);
		}

		switch(c)
		{
			case 0xEEE0: // normal
			{
				new_textspan(false, false, false, 0, token);
				break;
			}
			case 0xEEE1: // superscript
			{
				new_textspan(false, false, true, supbl, token);
				break;
			}
			case 0xEEE2: // subscript
			{
				new_textspan(false, false, true, subpl, token);
				break;
			}
			case 0xEEE4: // italic
			{
				new_textspan(false, true, false, 0, token);
				break;
			}
			case 0xEEE5: // italic superscript
			{
				new_textspan(false, true, true, supbl, token);
				break;
			}
			case 0xEEE6: // italic subscript
			{
				new_textspan(false, true, true, subpl, token);
				break;
			}
			case 0xEEE8: // bold normal
			{
				new_textspan(true, false, false, 0, token);
				break;
			}
			case 0xEEE9: // bold superscript
			{
				new_textspan(true, false, true, supbl, token);
				break;
			}
			case 0xEEEA: // bold subscript
			{
				new_textspan(true, false, true, subpl, token);
				break;
			}
			case 0xEEEC: // bold italic
			{
				new_textspan(true, true, false, 0, token);
				break;
			}
			case 0xEEED: // bold italic superscript
			{
				new_textspan(true, true, true, supbl, token);
				break;
			}
			case 0xEEEE: // bold italic subscript
			{
				new_textspan(true, true, true, subpl, token);
				break;
			}
			default:
			{
				prevchar = nextchar;
				continue;
			}
		}
		
		prevchar = nextchar;
	}
	while(search !== -1);
}

function toggle_calcdiv()
{
	var ele1 = document.getElementById("toggle_calcdiv_symbol");
	var ele2 = document.getElementById("calcdiv");
	var ele3 = document.getElementById("toggle_calcdiv_checkbox");
	
	if(!ele3.checked)
	{
		ele3.checked = true;
		ele2.style.display = "block";
		ele1.innerHTML = "\u25B2";
	}
	else
	{
		ele3.checked = false;
		ele2.style.display = "none";
		ele1.innerHTML = "\u25BC";
	}
}

function replace_minus_signs(str)
{
	return str.replace(/\u002D/g, "\u2212");
}

var decadetable = [
	10, 11, 12, 13, 14, 15, 16, 17, 18, 20,
	22, 24, 26, 28, 30, 32, 34, 36, 38, 40,
	45, 50, 55, 60, 65, 70, 75, 80, 85, 90
];

function expval(v, e, minexp, maxexp)
{
	if(newdecadescale)
		return decadetable[v % 30] * Math.pow(10, Math.floor(v / 30) - 10 + decadeshift);
	else
		return Math.pow(e, minexp + v * (maxexp - minexp) / 240 + decadeshift);
}

var expparams = [
	/*  0 */ [0, 0, 0],
	/*  1 */ [10, -1, 5],
	/*  2 */ [10, -7, -1],
	/*  3 */ [10, -9, -1],
	/*  4 */ [10, 0, 6],
	/*  5 */ [10, -9, -1],
	/*  6 */ [10, -5, -1],
	/*  7 */ [10, -9, -1],
	/*  8 */ [0, 0, 0],
	/*  9 */ [10, -9, -1],
	/* 10 */ [10, -9, -1],
	/* 11 */ [10, -9, 1],
	/* 12 */ [10, -5, 5],
	/* 13 */ [10, -11, -1],
	/* 14 */ [10, -12, -4],
	/* 15 */ [10, -12, -4]
];

function expval_wrap(v, index)
{
	return expval(v, expparams[index][0], expparams[index][1], expparams[index][2]);
}

function slider_input(index, noupdate, nocalculate)
{
	var slider = document.getElementById("slider" + index);
	var label = document.getElementById("value" + index);
	
	var rawval = (!slider ? undefined : Number(slider.value)); // [0,240]
	var val = 0;
	var valstr = "";
	var unitstr = "";
	
	switch(index)
	{
		case 3: // [S]_0
		{
			if(!nocalculate) S_0 = expval_wrap(rawval, index);
			val = S_0;
			unitstr = "mol\xA0l<sup>\u22121</sup>";
			break;
		}
		case 5: // [E]
		{
			if(!nocalculate) E_0 = expval_wrap(rawval, index);
			val = E_0;
			unitstr = "mol\xA0l<sup>\u22121</sup>";
			break;
		}
		case 7: // K_D
		{
			if(!nocalculate) K_D = expval_wrap(rawval, index);
			val = K_D;
			unitstr = "mol\xA0l<sup>\u22121</sup>";
			break;
		}
		case 8: // delta_G
		{
			val = 8.31446261815324 * 298.15 * Math.log(K_D) * 1e-3;
			unitstr = "kJ\xA0mol<sup>\u22121</sup>";
			break;
		}
		case 9: // K_D 2
		{
			if(!nocalculate) K_D2 = expval_wrap(rawval, index);
			val = K_D2;
			unitstr = "mol\xA0l<sup>\u22121</sup>";
			break;
		}
		case 10: // [E]
		{
			if(!nocalculate) Q_0 = expval_wrap(rawval, index);
			val = Q_0;
			unitstr = "mol\xA0l<sup>\u22121</sup>";
			break;
		}
		case 16: // K_A (association constant)
		{
			val = 1.0 / K_D;
			unitstr = "l\xA0mol<sup>\u22121</sup>";
			break;
		}
		case 17: // delta_G 2
		{
			val = 8.31446261815324 * 298.15 * Math.log(K_D2) * 1e-3;
			unitstr = "kJ\xA0mol<sup>\u22121</sup>";
			break;
		}
		case 18: // K_A 2 (association constant)
		{
			val = 1.0 / K_D2;
			unitstr = "l\xA0mol<sup>\u22121</sup>";
			break;
		}
	}
	
	var abs = Math.abs(val);
	var magnitude = Math.floor(Math.log10(abs));
	var sign = (val < 0) ? "\u2212" : "";
	
	if(magnitude === 0 || index === 17 || index === 8)
	{
		valstr = sign + abs.toFixed(value_digits);
	}
	else
	{
		valstr = sign + (abs / Math.pow(10, magnitude)).toFixed(value_digits) + "\xA0\u22C5\xA0" + "10<sup>" + magnitude + "</sup>";
	}
	
	label.innerHTML = replace_minus_signs(valstr) + "\xA0" + unitstr;
	
	if(index === 7)
	{
		slider_input(8, true);
		slider_input(16, true);
	}
	
	if(index === 9)
	{
		slider_input(17, true);
		slider_input(18, true);
	}
	
	if(!noupdate) update();
}

function radio_input(index, noupdate)
{
	switch(index)
	{
		case 1: // absolute y-axis
		{
			scale_absolute = 1;
			break;
		}
		case 2: // relative y-axis
		{
			scale_absolute = 0;
			break;
		}
		case 3: // x-axis type 1
		{
			xscale_alternative = 0;
			
			if(appmode === appmode_ligands)
			{
				document.getElementById("fixval5").disabled = true;
				document.getElementById("fixval10").disabled = false;
			}
			
			break;
		}
		case 4: // x-axis type 2
		{
			xscale_alternative = 1;
			
			if(appmode === appmode_ligands)
			{
				document.getElementById("fixval5").disabled = false;
				document.getElementById("fixval10").disabled = true;
			}
			
			break;
		}
		case 11: // x-axis type 3
		{
			xscale_alternative = 2;
			
			if(appmode === appmode_receptors)
			{
				document.getElementById("fixval5").disabled = false;
				document.getElementById("fixval10").disabled = true;
			}
			
			break;
		}
		case 5: // absolute logarithmic y-axis
		{
			scale_absolute = 1;
			break;
		}
		case 6: // absolute linear y-axis
		{
			scale_absolute = 0;
			break;
		}
		case 7: // specificity axis
		{
			scale_absolute = 2;
			break;
		}
		default:
		{
			break;
		}
	}
	
	if(!noupdate) update();
}

function remove_from_parent(ele)
{
	if(!ele) return null;
	return ele.parentNode.removeChild(ele);
}

function delete_all_children(ele)
{
	if(!ele) return;
	while(ele.hasChildNodes()) ele.removeChild(ele.firstChild);
}

var resizing = false;

function mouse_moved()
{
	var figure_div = document.getElementById("figure_div");
	var fwidth = extract_px(figure_div.style.width);
	var fheight = extract_px(figure_div.style.height);
	
	if(!(fwidth === fwidth_old && fheight === fheight_old))
	{
		fwidth_old = fwidth;
		fheight_old = fheight;
		resizing = true;
		update(false);
	}
}

function mouse_released()
{
	if(resizing)
	{
		resizing = false;
		update();
	}
}

function data_changed(do_update)
{
	var databox = document.getElementById("databox");
	var lines = [databox.value];
	var n, p, num1, num2;
	var line, part1, part2;
	var datastr = "";
	
	datapoints = [];
	
	while(1)
	{
		n = lines.length - 1;
		p = lines[n].indexOf("\n");
		
		if(p === -1) break;
		
		part1 = lines[n].substring(0, p);
		part2 = lines[n].substring(p + 1);
		
		lines[n] = part1;
		lines[n + 1] = part2;
	}

	for(n = 0; n < lines.length; n++)
	{
		line = lines[n].trim();
		p = line.lastIndexOf("\x20");
		
		if(p === -1)
		{
			p = line.lastIndexOf("\x09");
			if(p === -1) continue;
		}
		
		part1 = line.substring(0, p + 1).trim();
		part2 = line.substring(p + 1).trim();
		
		num1 = Number(part1);
		num2 = Number(part2);
		
		if(part1.length === 0 || part2.length === 0 || Number.isNaN(num1) || Number.isNaN(num2)) continue;
		
		datapoints[datapoints.length] = {x: num1, y: num2};
		
		// datastr += num1.toExponential() + "\x09" + num2.toExponential() + "\n";
	}
	
	// databox.value = datastr;
	
	delete_all_children(document.getElementById("svg_datapointgroup"));
	
	if(do_update === undefined || do_update)
	{
		// calculate_lsf();
		update(false);
	}
}

function extract_px(str)
{
	return Number(str.substring(0, str.search("px")));
}

function zeropad(number, digits)
{
	var str = number.toFixed(0);
	return "0".repeat(str.length < digits ? digits - str.length : 0) + str;
}

function export_svg()
{
	var figure_div = document.getElementById("figure_div");
	
	var date = new Date(Date.now());
	var filename = "graphic_" + date.getFullYear() + "-" + zeropad(date.getMonth() + 1, 2) + "-" + zeropad(date.getDate(), 2) + "_" + 
			zeropad(date.getHours(), 2) + "-" + zeropad(date.getMinutes(), 2) + "-" + zeropad(date.getSeconds(), 2) + ".svg";
	
	var str = figure_div.innerHTML.replace("width=\"96%\" ", "").replace("height=\"100%\" ", "").replace(/\t/g, "").replace(/\x3C\x2Fpolyline\x3E\x3Cpolyline/g, "</polyline>\n\t<polyline");
	var blob = new File(["<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"no\" ?>\n", str], filename, {type: "application/svg+xml"});
	var url = URL.createObjectURL(blob);
	window.open(url, "_blank");
}
