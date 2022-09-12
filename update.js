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

function update(recalculate)
{
	var figure = document.getElementById("figure_svg");
	var figure_piegroup = document.getElementById("svg_piegroup");
	var figure_curvegroup = document.getElementById("svg_curvegroup");
	var figure_datapointgroup = document.getElementById("svg_datapointgroup");
	var figure_overlaygroup = document.getElementById("svg_overlaygroup");
	var figure_linegroup = document.getElementById("svg_linegroup");
	var figure_legendgroup = document.getElementById("svg_legendgroup");
	var figure_scribblegroup = document.getElementById("svg_scribblegroup");
	
	var figure_div = document.getElementById("figure_div");
	var cwidth = extract_px(figure_div.style.width) * 0.96;
	var cheight = extract_px(figure_div.style.height);
	var cminwidth = extract_px(figure_div.style.minWidth) * 0.96;
	var cminheight = extract_px(figure_div.style.minHeight);

	if(cwidth < cminwidth)
	{
		figure_div.style.width = figure_div.style.minWidth;
		cwidth = cminwidth;
	}
	
	if(cheight < cminheight)
	{
		figure_div.style.height = figure_div.style.minHeight;
		cheight = cminheight;
	}

	figure.setAttribute("viewBox", "0 0 " + cwidth + " " + cheight);
	
	if(recalculate === undefined) recalculate = true;
	
	var ca_left = 60;
	var ca_right = 60;
	var ca_top = 20;
	var ca_bottom = 40;
	var ca_width = cwidth - ca_left - ca_right;
	var ca_height = cheight - ca_top - ca_bottom;
	
	var c = 0;
	var v = 0;
	var h = 0;
	var t = 0;
	var S = S_0;
	var S_eq = 0;
	var S_eff = S;
	
	var cd = undefined;
	var cm = undefined;
	var masses = undefined;
	
	var ele_id = "", ele2_id = "";
	
	var colour1 = "rgb(255,0,0)";
	var colour2 = "rgb(0,0,255)";
	var colour3 = "rgb(0,128,0)";
	var colour4 = "rgb(128,0,128)";
	var colour5 = "rgb(170,85,0)";
	var colour6 = "rgb(255,128,0)";
	var colour7 = "rgb(128,128,128)";
	
	if(recalculate) switch(appmode)
	{
		case appmode_ligand:
		{
			curves[0] = [];
			curves[1] = [];
			curves[2] = [];
			curves[3] = [];
			curves[4] = [];
			
			labels[0] = "\uEEE8[PL]";
			labels[1] = "\uEEE8[P]";
			labels[3] = "\uEEECc\uEEEAP";
			labels[4] = "\uEEECK\uEEEAD";
			
			legends[0] = 0;
			legends[1] = 1;
			legends[2] = null;
			legends[3] = null;
			legends[4] = null;
			
			colours[0] = colour1;
			colours[1] = colour2;
			colours[2] = colour5;
			colours[3] = colour6;
			colours[4] = colour6;

			if(xscale_alternative) // free ligand
			{
				for(h = 0; h <= 20; h += h < 1 ? 1/64 : 1/16)
				{
					c = (h === 0 ? 1e-10 : h) * K_D;
					cd = calculate_ligand_free(E_0, c, K_D).d;
					
					curves[0].push({
						x: c,
						y: (scale_absolute ? Math.log10(cd[0]) : (cd[0] / E_0))
					});
					
					curves[1].push({
						x: c,
						y: (scale_absolute ? Math.log10(cd[1]) : (cd[1] / E_0))
					});
				}
			}
			else // total ligand
			{
				for(h = decadeshift - 10; h <= decadeshift; h += 1/32)
				{
					c = Math.pow(10, h);
					cd = calculate_ligand_total(E_0, c, K_D).d;
					
					curves[0].push({
						x: h,
						y: (scale_absolute ? Math.log10(cd[0]) : (cd[0] / E_0))
					});
					
					curves[1].push({
						x: h,
						y: (scale_absolute ? Math.log10(cd[1]) : (cd[1] / E_0))
					});
				}
			}
			
			if(scale_absolute)
			{
				ymin = decadeshift - 12;
				ymax = decadeshift;
				yaxistype = axistype_log;
			}
			else
			{
				ymin = 0;
				ymax = 1;
				yaxistype = axistype_lin;
			}
			
			var datapoint = calculate_ligand_total(E_0, S_0, K_D);
			cd = datapoint.d;
			cm = datapoint.m;
			
			var m1 = Number(document.getElementById("mass1").value);
			var m2 = Number(document.getElementById("mass2").value);
			if(!(m1 > 0)) m1 = NaN;
			if(!(m2 > 0)) m2 = NaN;
			masses = [m1 + m2, m1];
			
			pie[0] = { value: (cd[0] / E_0), colour: colour1 };
			pie[1] = { value: (cd[1] / E_0), colour: colour2 };
			
			if(xscale_alternative)
			{
				xmin = 0;
				xmax = 20 * K_D;
				xaxistype = axistype_lin;
				
				if(cd[2] >= 0 && cd[2] <= 1.0001 * xmax)
				{
					curves[2][0] = {x: cd[2], y: ymin};
					curves[2][1] = {x: cd[2], y: ymax + 6 * (ymax - ymin) / ca_height};
				}
				else
				{
					curves[2] = undefined;
				}
				
				labels[2] = "\uEEE8[L]";
				
				if(E_0 <= 1.0001 * xmax)
				{
					curves[3][0] = {x: E_0, y: ymin};
					curves[3][1] = {x: E_0, y: ymax + 6 * (ymax - ymin) / ca_height};
				}
				else
				{
					curves[3] = undefined;
				}
				
				if(K_D <= 1.0001 * xmax)
				{
					curves[4][0] = {x: K_D, y: ymin};
					curves[4][1] = {x: K_D, y: ymax + 6 * (ymax - ymin) / ca_height};
				}
				else
				{
					curves[4] = undefined;
				}
			}
			else
			{
				xmin = decadeshift - 10;
				xmax = decadeshift;
				xaxistype = axistype_log;
				
				if(Math.log10(S_0) >= xmin && Math.log10(S_0) <= xmax)
				{
					curves[2][0] = {x: Math.log10(S_0), y: ymin};
					curves[2][1] = {x: Math.log10(S_0), y: ymax + 6 * (ymax - ymin) / ca_height};
				}
				else
				{
					curves[2] = undefined;
				}
				
				labels[2] = "\uEEECc\uEEEAL";
				
				if(Math.log10(E_0) >= xmin && Math.log10(E_0) <= xmax)
				{
					curves[3][0] = {x: Math.log10(E_0), y: ymin};
					curves[3][1] = {x: Math.log10(E_0), y: ymax + 6 * (ymax - ymin) / ca_height};
				}
				else
				{
					curves[3] = undefined;
				}
				
				if(Math.log10(K_D) >= xmin && Math.log10(K_D) <= xmax)
				{
					curves[4][0] = {x: Math.log10(K_D), y: ymin};
					curves[4][1] = {x: Math.log10(K_D), y: ymax + 6 * (ymax - ymin) / ca_height};
				}
				else
				{
					curves[4] = undefined;
				}
			}
			
			break;
		}
		case appmode_homodimer:
		{
			curves[0] = [];
			curves[1] = [];
			curves[2] = [];
			curves[3] = [];
			
			labels[0] = "\uEEE8[P\uEEEA2\uEEE8]";
			labels[1] = "\uEEE8[P]";
			labels[3] = "\uEEECK\uEEEAD";
			
			legends[0] = 0;
			legends[1] = 1;
			legends[2] = null;
			legends[3] = null;
			
			colours[0] = colour1;
			colours[1] = colour2;
			colours[2] = colour5;
			colours[3] = colour6;
			
			if(xscale_alternative) // free protein
			{
				for(h = 0; h <= 10; h += h < 0.5 ? 1/128 : 1/32)
				{
					c = (h === 0 ? 1e-10 : h) * K_D;
					cd = calculate_homodimer_free(c, K_D).d;
					
					curves[0].push({
						x: c,
						y: (scale_absolute ? Math.log10(cd[0]) : 2 * cd[0] / (2 * cd[0] + cd[1]))
					});
					
					curves[1].push({
						x: c,
						y: (scale_absolute ? Math.log10(cd[1]) : cd[1] / (2 * cd[0] + cd[1]))
					});
				}
				
				labels[2] = "\uEEE8[P]";
			}
			else // total protein
			{
				for(h = decadeshift - 10; h <= decadeshift; h += 1/32)
				{
					c = Math.pow(10, h);
					cd = calculate_homodimer(c, K_D).d;
					
					curves[0].push({
						x: h,
						y: (scale_absolute ? Math.log10(cd[0]) : 2 * cd[0] / (2 * cd[0] + cd[1]))
					});
					
					curves[1].push({
						x: h,
						y: (scale_absolute ? Math.log10(cd[1]) : cd[1] / (2 * cd[0] + cd[1]))
					});
				}
				
				labels[2] = "\uEEECc\uEEEAP";
			}
			
			if(scale_absolute)
			{
				ymin = decadeshift - 12;
				ymax = decadeshift;
				yaxistype = axistype_log;
			}
			else
			{
				ymin = 0;
				ymax = 1;
				yaxistype = axistype_lin;
			}
			
			var datapoint = calculate_homodimer(E_0, K_D);
			cd = datapoint.d;
			cm = datapoint.m;
			
			var m1 = Number(document.getElementById("mass1").value);
			if(!(m1 > 0)) m1 = NaN;
			masses = [m1 * 2, m1];
			
			pie[0] = {value: 2 * cd[0] / E_0, colour: colour1};
			pie[1] = {value: cd[1] / E_0, colour: colour2};
			
			if(xscale_alternative)
			{
				xmin = 0;
				xmax = 10 * K_D;
				xaxistype = axistype_lin;
				
				if(cd[1] <= 1.0001 * xmax)
				{
					curves[2][0] = {x: cd[1], y: ymin};
					curves[2][1] = {x: cd[1], y: ymax + 6 * (ymax - ymin) / ca_height};
				}
				else
				{
					curves[2] = undefined;
				}
				
				if(K_D <= 1.0001 * xmax)
				{
					curves[3][0] = {x: K_D, y: ymin};
					curves[3][1] = {x: K_D, y: ymax + 6 * (ymax - ymin) / ca_height};
				}
				else
				{
					curves[3] = undefined;
				}
			}
			else
			{
				xmin = decadeshift - 10;
				xmax = decadeshift;
				xaxistype = axistype_log;
				
				if(Math.log10(E_0) >= xmin && Math.log10(E_0) <= xmax)
				{ 
					curves[2][0] = {x: Math.log10(E_0), y: ymin};
					curves[2][1] = {x: Math.log10(E_0), y: ymax + 6 * (ymax - ymin) / ca_height};
				}
				else
				{
					curves[2] = undefined;
				}
				
				if(Math.log10(K_D) >= xmin && Math.log10(K_D) <= xmax)
				{ 
					curves[3][0] = {x: Math.log10(K_D), y: ymin};
					curves[3][1] = {x: Math.log10(K_D), y: ymax + 6 * (ymax - ymin) / ca_height};
				}
				else
				{
					curves[3] = undefined;
				}
			}
			
			break;
		}
		case appmode_ligands:
		{
			curves[0] = [];
			curves[1] = [];
			curves[2] = [];
			curves[3] = [];
			curves[4] = [];
			curves[5] = [];
			curves[6] = [];
			curves[7] = [];
			curves[8] = [];
			
			labels[0] = "\uEEE8[PL]";
			labels[1] = "\uEEE8[PL\u2032]";
			labels[2] = "\uEEE8[P]";
			labels[3] = "";
			labels[4] = "";
			labels[5] = "";
			labels[6] = "\uEEECK\uEEEAD";
			labels[7] = "\uEEECK\uEEEAD\uEEE8\u2032";
			labels[8] = "\uEEECc\uEEEAL\u2032\uEEE8";
			
			legends[0] = 0;
			legends[1] = 1;
			legends[2] = 2;
			legends[3] = null;
			legends[4] = null;
			legends[5] = null;
			legends[6] = null;
			legends[7] = null;
			legends[8] = null;
			
			colours[0] = colour1;
			colours[1] = colour2;
			colours[2] = colour3;
			colours[3] = colour4;
			colours[4] = colour6;
			colours[5] = colour6;
			colours[6] = colour6;
			colours[7] = colour6;
			colours[8] = colour6;
			
			if(xscale_alternative)
				colours[4] = colour5;
			else
				colours[5] = colour5;
			
			for(h = decadeshift - 10; h <= decadeshift; h += 1/32)
			{
				c = Math.pow(10, h);
				
				if(xscale_alternative) // x-axis is ligand concentration
					cd = calculate_ligands(E_0, c, S_0, K_D, K_D2).d;
				else // x-axis is protein concentration
					cd = calculate_ligands(c, Q_0, S_0, K_D, K_D2).d;
				
				curves[0].push({
					x: h,
					y: (scale_absolute ? Math.log10(cd[0]) : cd[0] / (cd[0] + cd[1] + cd[2]))
				});
				
				curves[1].push({
					x: h,
					y: (scale_absolute ? Math.log10(cd[1]) : cd[1] / (cd[0] + cd[1] + cd[2]))
				});
				
				curves[2].push({
					x: h,
					y: (scale_absolute ? Math.log10(cd[2]) : cd[2] / (cd[0] + cd[1] + cd[2]))
				});
			}
			
			xmin = decadeshift - 10;
			xmax = decadeshift;
			xaxistype = axistype_log;
			
			if(scale_absolute)
			{
				ymin = decadeshift - (extmode ? 12 : 10);
				ymax = decadeshift;
				yaxistype = axistype_log;
			}
			else
			{
				ymin = 0;
				ymax = 1;
				yaxistype = axistype_lin;
			}
			
			if(Math.log10(Q_0) >= xmin && Math.log10(Q_0) <= xmax)
			{
				curves[4][0] = {x: Math.log10(Q_0), y: ymin};
				curves[4][1] = {x: Math.log10(Q_0), y: ymax + 6 * (ymax - ymin) / ca_height};
			}
			else
			{
				curves[4] = undefined;
			}
			
			labels[4] = "\uEEECc\uEEEAL";
			
			if(Math.log10(E_0) >= xmin && Math.log10(E_0) <= xmax)
			{
				curves[5][0] = {x: Math.log10(E_0), y: ymin};
				curves[5][1] = {x: Math.log10(E_0), y: ymax + 6 * (ymax - ymin) / ca_height};
			}
			else
			{
				curves[5] = undefined;
			}
			
			labels[5] = "\uEEECc\uEEEAP";
			
			if(Math.log10(K_D) >= xmin && Math.log10(K_D) <= xmax)
			{
				curves[6][0] = {x: Math.log10(K_D), y: ymin};
				curves[6][1] = {x: Math.log10(K_D), y: ymax + 6 * (ymax - ymin) / ca_height};
			}
			else
			{
				curves[6] = undefined;
			}
			
			if(Math.log10(K_D2) >= xmin && Math.log10(K_D2) <= xmax)
			{
				curves[7][0] = {x: Math.log10(K_D2), y: ymin};
				curves[7][1] = {x: Math.log10(K_D2), y: ymax + 6 * (ymax - ymin) / ca_height};
			}
			else
			{
				curves[7] = undefined;
			}
			
			if(Math.log10(S_0) >= xmin && Math.log10(S_0) <= xmax)
			{
				curves[8][0] = {x: Math.log10(S_0), y: ymin};
				curves[8][1] = {x: Math.log10(S_0), y: ymax + 6 * (ymax - ymin) / ca_height};
			}
			else
			{
				curves[8] = undefined;
			}
			
			var datapoint = calculate_ligands(E_0, Q_0, S_0, K_D, K_D2);
			cd = datapoint.d;
			cm = datapoint.m;
			
			pie[0] = {value: cd[0] / (cd[0] + cd[1] + cd[2]), colour: colour1};
			pie[1] = {value: cd[1] / (cd[0] + cd[1] + cd[2]), colour: colour2};
			pie[2] = {value: cd[2] / (cd[0] + cd[1] + cd[2]), colour: colour3};
			
			var m1 = Number(document.getElementById("mass1").value);
			var m2 = Number(document.getElementById("mass2").value);
			var m3 = Number(document.getElementById("mass3").value);
			if(!(m1 > 0)) m1 = NaN;
			if(!(m2 > 0)) m2 = NaN;
			if(!(m3 > 0)) m3 = NaN;
			masses = [m1 + m2, m1 + m3, m1];
			
			break;
		}
		case appmode_receptors:
		{
			curves[0] = [];
			curves[1] = [];
			curves[2] = [];
			curves[3] = [];
			curves[4] = [];
			curves[5] = [];
			curves[6] = [];
			curves[7] = [];
			curves[8] = [];
			
			labels[0] = "\uEEE8[PL]";
			labels[1] = "\uEEE8[P\u2032L]";
			labels[2] = "\uEEE8[P]";
			labels[3] = "\uEEE8[P\u2032]";
			labels[4] = "\uEEECc\uEEEAP";
			labels[5] = "\uEEECc\uEEEAL";
			labels[6] = "\uEEECK\uEEEAD";
			labels[7] = "\uEEECK\uEEEAD\uEEE8\u2032";
			labels[8] = "\uEEECc\uEEEAP\u2032\uEEE8";
			
			legends[0] = 0;
			legends[1] = 1;
			legends[2] = 2;
			legends[3] = 3;
			legends[4] = null;
			legends[5] = null;
			legends[6] = null;
			legends[7] = null;
			legends[8] = null;
			
			colours[0] = colour1;
			colours[1] = colour2;
			colours[2] = colour3;
			colours[3] = colour4;
			colours[4] = colour6;
			colours[5] = colour5;
			colours[6] = colour6;
			colours[7] = colour6;
			colours[8] = colour6;
			
			if(scale_absolute === 2)
			{
				curves[1] = curves[2] = curves[3] = undefined;
				legends[1] = legends[2] = legends[3] = 0;
				labels[1] = labels[2] = labels[3] = undefined;
				colours[0] = colour7;
				labels[0] = "\uEEEC\u03B1\uEEEAs";
			}
			
			if(xscale_alternative) // x-axis is free ligand concentration
			{
				xmin = 0;
				xmax = 20 * Math.min(K_D, K_D2);
				xaxistype = axistype_lin;
				
				labels[5] = "\uEEE8[L]";
				
				for(h = 0; h <= 20; h += h < 1 ? 1/64 : 1/16)
				{
					c = (h === 0 ? 1e-10 : h) * xmax / 20;
					cd = calculate_ligands_free(c, Q_0, S_0, K_D, K_D2).d;
					
					if(scale_absolute === 2)
					{
						curves[0].push({
							x: c,
							y: cd[5] // specificity
						});
					}
					else
					{
						curves[0].push({
							x: c,
							y: (scale_absolute ? Math.log10(cd[0]) : cd[0])
						});
						
						curves[1].push({
							x: c,
							y: (scale_absolute ? Math.log10(cd[1]) : cd[1])
						});
						
						curves[2].push({
							x: c,
							y: (scale_absolute ? Math.log10(cd[2]) : cd[2])
						});
						
						curves[3].push({
							x: c,
							y: (scale_absolute ? Math.log10(cd[3]) : cd[3])
						});
					}
				}
				
				var wmax = function(a,b) { return Math.max(a, b.y); };
				
				if(scale_absolute === 1)
				{
					ymin = decadeshift - (extmode ? 12 : 10);
					ymax = decadeshift;
					yaxistype = axistype_log;
				}
				else if(scale_absolute === 2)
				{
					ymin = 0;
					ymax = curves[0].reduce(wmax, 0);
					yaxistype = axistype_lin;
				}
				else
				{
					ymin = 0;
					ymax = Math.max(curves[0].reduce(wmax, 0), curves[1].reduce(wmax, 0));
					yaxistype = axistype_lin;
				}
				
				if(Q_0 <= 1.0001 * xmax)
				{
					curves[4][0] = {x: Q_0, y: ymin};
					curves[4][1] = {x: Q_0, y: ymax + 6 * (ymax - ymin) / ca_height};
				}
				else
				{
					curves[4] = undefined;
				}
				
				if(K_D <= 1.0001 * xmax)
				{
					curves[6][0] = {x: K_D, y: ymin};
					curves[6][1] = {x: K_D, y: ymax + 6 * (ymax - ymin) / ca_height};
				}
				else
				{
					curves[6] = undefined;
				}
				
				if(K_D2 <= 1.0001 * xmax)
				{
					curves[7][0] = {x: K_D2, y: ymin};
					curves[7][1] = {x: K_D2, y: ymax + 6 * (ymax - ymin) / ca_height};
				}
				else
				{
					curves[7] = undefined;
				}
				
				if(S_0 <= 1.0001 * xmax)
				{
					curves[8][0] = {x: S_0, y: ymin};
					curves[8][1] = {x: S_0, y: ymax + 6 * (ymax - ymin) / ca_height};
				}
				else
				{
					curves[8] = undefined;
				}
			}
			else
			{
				xmin = decadeshift - 10;
				xmax = decadeshift;
				xaxistype = axistype_log;
				
				labels[5] = "\uEEECc\uEEEAL";
				
				for(h = decadeshift - 10; h <= decadeshift; h += 1/32)
				{
					c = Math.pow(10, h);
					cd = calculate_ligands(c, Q_0, S_0, K_D, K_D2).d;
					
					if(scale_absolute === 2)
					{
						curves[0].push({
							x: h,
							y: cd[5] // specificity
						});
					}
					else
					{
						curves[0].push({
							x: h,
							y: (scale_absolute ? Math.log10(cd[0]) : cd[0])
						});
						
						curves[1].push({
							x: h,
							y: (scale_absolute ? Math.log10(cd[1]) : cd[1])
						});
						
						curves[2].push({
							x: h,
							y: (scale_absolute ? Math.log10(cd[2]) : cd[2])
						});
						
						curves[3].push({
							x: h,
							y: (scale_absolute ? Math.log10(cd[3]) : cd[3])
						});
					}
				}
				
				var wmax = function(a,b) { return Math.max(a, b.y); };
				
				if(scale_absolute === 1)
				{
					ymin = decadeshift - (extmode ? 12 : 10);
					ymax = decadeshift;
					yaxistype = axistype_log;
				}
				else if(scale_absolute === 2)
				{
					ymin = 0;
					ymax = curves[0].reduce(wmax, 0);
					yaxistype = axistype_lin;
				}
				else
				{
					ymin = 0;
					ymax = Math.max(curves[0].reduce(wmax, 0), curves[1].reduce(wmax, 0), curves[2].reduce(wmax, 0), curves[3].reduce(wmax, 0));
					yaxistype = axistype_lin;
				}
				
				if(Math.log10(Q_0) >= xmin && Math.log10(Q_0) <= xmax)
				{
					curves[4][0] = {x: Math.log10(Q_0), y: ymin};
					curves[4][1] = {x: Math.log10(Q_0), y: ymax + 6 * (ymax - ymin) / ca_height};
				}
				else
				{
					curves[4] = undefined;
				}
				
				if(Math.log10(K_D) >= xmin && Math.log10(K_D) <= xmax)
				{
					curves[6][0] = {x: Math.log10(K_D), y: ymin};
					curves[6][1] = {x: Math.log10(K_D), y: ymax + 6 * (ymax - ymin) / ca_height};
				}
				else
				{
					curves[6] = undefined;
				}
				
				if(Math.log10(K_D2) >= xmin && Math.log10(K_D2) <= xmax)
				{
					curves[7][0] = {x: Math.log10(K_D2), y: ymin};
					curves[7][1] = {x: Math.log10(K_D2), y: ymax + 6 * (ymax - ymin) / ca_height};
				}
				else
				{
					curves[7] = undefined;
				}
				
				if(Math.log10(S_0) >= xmin && Math.log10(S_0) <= xmax)
				{
					curves[8][0] = {x: Math.log10(S_0), y: ymin};
					curves[8][1] = {x: Math.log10(S_0), y: ymax + 6 * (ymax - ymin) / ca_height};
				}
				else
				{
					curves[8] = undefined;
				}
			}
			
			var datapoint;
			
			datapoint = calculate_ligands(E_0, Q_0, S_0, K_D, K_D2);
			
			cd = datapoint.d;
			cm = datapoint.m;
			
			pie[0] = {value: cd[0] / (cd[0] + cd[1] + cd[2] + cd[3]), colour: colour1};
			pie[1] = {value: cd[1] / (cd[0] + cd[1] + cd[2] + cd[3]), colour: colour2};
			pie[2] = {value: cd[2] / (cd[0] + cd[1] + cd[2] + cd[3]), colour: colour3};
			pie[3] = {value: cd[3] / (cd[0] + cd[1] + cd[2] + cd[3]), colour: colour4};
			
			var m1 = Number(document.getElementById("mass1").value);
			var m2 = Number(document.getElementById("mass2").value);
			var m3 = Number(document.getElementById("mass3").value);
			if(!(m1 > 0)) m1 = NaN;
			if(!(m2 > 0)) m2 = NaN;
			if(!(m3 > 0)) m3 = NaN;
			masses = [m1 + m2, m1 + m3, m2, m3];
			
			if(xscale_alternative)
			{
				if(cd[4] >= 0 && cd[4] <= 1.0001 * xmax)
				{
					curves[5][0] = {x: cd[4], y: ymin};
					curves[5][1] = {x: cd[4], y: ymax + 6 * (ymax - ymin) / ca_height};
				}
				else
				{
					curves[5] = undefined;
				}
			}
			else
			{
				if(Math.log10(E_0) >= xmin && Math.log10(E_0) <= xmax)
				{
					curves[5][0] = {x: Math.log10(E_0), y: ymin};
					curves[5][1] = {x: Math.log10(E_0), y: ymax + 6 * (ymax - ymin) / ca_height};
				}
				else
				{
					curves[5] = undefined;
				}
			}
			
			break;
		}
	}
	
	// #################################################################################################################
	
	if(cd !== undefined && cm !== undefined)
	{
		var i;
		var csum = 0;
		var concstr, epos;
		
		for(i = 0; i < datalabels[appmode].length; i++)
		{
			csum += cd[i] * cm[i];
		}
		
		for(i = 0; i < datalabels[appmode].length; i++)
		{
			var ele1 = document.getElementById("data" + i + "l");
			var ele2 = document.getElementById("data" + i + "a");
			var ele3 = document.getElementById("data" + i + "b");
			var ele4 = document.getElementById("data" + i + "c");
			
			concstr = cd[i].toExponential(value_digits);
			epos = concstr.search("e");
			
			if(epos !== -1)
			{
				if(Number(concstr.substring(epos + 1)) === 0)
				{
					ele2.innerHTML = replace_minus_signs(concstr.substring(0, epos));
				}
				else
				{
					ele2.innerHTML = replace_minus_signs(concstr.substring(0, epos) +
						"\xA0\u22C5\xA010<sup>" +
						Number(concstr.substring(epos + 1)).toString() +
						"</sup>");
				}
			}
			else
			{
				ele2.innerHTML = replace_minus_signs(concstr);
			}
			
			concstr = (masses[i] * cd[i]).toExponential(value_digits);
			epos = concstr.search("e");
			
			if(!(masses[i] > 0))
			{
				ele3.innerHTML = "&mdash;";
			}
			else if(epos !== -1)
			{
				if(Number(concstr.substring(epos + 1)) === 0)
				{
					ele3.innerHTML = replace_minus_signs(concstr.substring(0, epos));
				}
				else
				{
					ele3.innerHTML = replace_minus_signs(concstr.substring(0, epos) +
						"\xA0\u22C5\xA010<sup>" +
						Number(concstr.substring(epos + 1)).toString() +
						"</sup>");
				}
			}
			else
			{
				ele3.innerHTML = replace_minus_signs(concstr);
			}
			
			ele4.innerHTML = replace_minus_signs((100 * cd[i] * cm[i] / csum).toFixed(1));
			
			if(appmode === appmode_receptors)
			{
				var csum2 = 0;
				var colourstr = "";
				
				if(i === 0 || i === 2)
				{
					csum2 = cd[0] + cd[2];
					colourstr = colour1;
				}
				else
				{
					csum2 = cd[1] + cd[3];
					colourstr = colour2;
				}
				
				ele4.innerHTML += " / <span style=\"color:" + colourstr + "\">" +
					replace_minus_signs((100 * cd[i] * cm[i] / csum2).toFixed(1)) + "</span>";
			}
			
			ele1.style.color = (i === 0 && scale_absolute === 2) ? colour1 : colours[i];
		}
		
		var specificity_label = document.getElementById("specificity_label");

		if(specificity_label)
		{
			var alpha_s = cd[0] / cd[1];
			
			concstr = alpha_s.toExponential(value_digits);
			epos = concstr.search("e");
			
			if(epos !== -1)
			{
				if(Number(concstr.substring(epos + 1)) === 0)
				{
					specificity_label.innerHTML = replace_minus_signs(concstr.substring(0, epos));
				}
				else
				{
					specificity_label.innerHTML = replace_minus_signs(concstr.substring(0, epos) +
						"\xA0\u22C5\xA010<sup>" +
						Number(concstr.substring(epos + 1)).toString() +
						"</sup>");
				}
			}
			else
			{
				specificity_label.innerHTML = replace_minus_signs(concstr);
			}
		}
	}
	
	function transform(point)
	{
		if(point)
		{
			var point2 = {
				x: ca_width * (point.x - xmin) / (xmax - xmin),
				y: ca_height * (point.y - ymax) / (ymin - ymax)
			};
			return point2;
		}
		else
		{
			return undefined;
		}
	}
	
	var p = {x:0, y:0};

	if(pie.length > 0)
	{
		var radius = (cwidth - ca_left - ca_width) * 0.5 - 2;
		var cx = ca_left + ca_width + radius + 2;
		var cy = ca_top + radius + 45 + 40 * legends.reduce(function(max, val) { return ((val > 0) && val > max) ? val : max; }, 0);
		var sigma_value = 0;
		
		for(h = 0; h < pie.length; h++)
		{
			var x1 = cx + radius * Math.cos(-(0.5 - 2 * sigma_value) * Math.PI);
			var y1 = cy + radius * Math.sin(-(0.5 - 2 * sigma_value) * Math.PI);
			var x2 = cx + radius * Math.cos(-(0.5 - 2 * (sigma_value + pie[h].value)) * Math.PI);
			var y2 = cy + radius * Math.sin(-(0.5 - 2 * (sigma_value + pie[h].value)) * Math.PI);
			
			ele_id = "svg_pie" + h;
			var svg_arc = document.getElementById(ele_id);
			
			if(pie[h].value === 0)
			{
				if(svg_arc) remove_from_parent(svg_arc);
			}
			else
			{
				if(!svg_arc)
				{
					svg_arc = document.createElementNS(svg_xmlns, "path");
					svg_arc.style.fill = pie[h].colour;
					svg_arc.style.stroke = "none";
					svg_arc.id = ele_id;
					figure_piegroup.appendChild(svg_arc);
				}
				
				if(pie[h].value < 0.99)
				{
					svg_arc.setAttribute("d", "M " + cx + " " + cy + " L " + x1 + " " + y1 + " A " + radius + " " + radius + " 0 " + (pie[h].value > 0.5 ? "1" : "0") + " 1 " + x2 + " " + y2); 
				}
				else
				{
					var x3 = cx + radius * Math.cos(-(1.3 - 2 * sigma_value) * Math.PI);
					var y3 = cy + radius * Math.sin(-(1.3 - 2 * sigma_value) * Math.PI);
					svg_arc.setAttribute("d", "M " + cx + " " + cy + " L " + x1 + " " + y1 + " A " + radius + " " + radius + " 0 1 1 " + x3 + " " + y3 + " A " + radius + " " + radius + " 0 0 1 " + x2 + " " + y2); 
				}
				
				if(extpiemode)
				{
					var ele2_id = "svg_piesectormark" + h;
					var ele3_id = "svg_piesectorlabel" + h;
					var svg_sm = document.getElementById(ele2_id);
					var svg_sl = document.getElementById(ele3_id);
					
					if(!svg_sm)
					{
						svg_sm = document.createElementNS(svg_xmlns, "path");
						svg_sm.style.fill = "none";
						svg_sm.style.stroke = pie[h].colour;
						svg_sm.style.strokeWidth = "3px";
						svg_sm.style.strokeLinejoin = "round";
						svg_sm.style.strokeLinecap = "butt";
						svg_sm.id = ele2_id;
						figure_scribblegroup.appendChild(svg_sm);
					}
					
					if(!svg_sl)
					{
						svg_sl = document.createElementNS(svg_xmlns, "text");
						svg_sl.style.fill = pie[h].colour;
						svg_sl.setAttributeNS(xml_xmlns, "xml:space", "preserve");
						svg_sl.style.textAnchor = "end";
						svg_sl.id = ele3_id;
						figure_scribblegroup.appendChild(svg_sl);
					}
					
					var tv;
					var centre = false;
					
					if(pie.length === 2 && h < pie.length - 1)
					{
						if(pie[0].value > 2/3)
						{
							tv = 0.5 - (pie[1].value / 2 + pie[0].value);
						}
						else if(pie[0].value < 1/6)
						{
							tv = sigma_value + pie[h].value / 2;
						}
						else
						{
							tv = sigma_value + pie[h].value / 2;
							centre = true;
						}
					}
					else if(pie.length === 4 && h === 2)
					{
						tv = 0.5 - sigma_value;
					}
					else
					{
						tv = sigma_value + pie[h].value / 2;
					}
					
					if(pie.length === 3 && h === 0) centre = true;
					if(pie.length === 4 && h === 0) centre = true;
					
					var x4 = cx + (centre ? radius * 2 / 3 : radius + 3) * Math.cos(-(0.5 - 2 * tv) * Math.PI);
					var y4 = cy + (centre ? radius * 2 / 3 : radius + 3) * Math.sin(-(0.5 - 2 * tv) * Math.PI);
					var x5 = cx + (centre ? radius * 2 / 3 : radius + 8) * Math.cos(-(0.5 - 2 * tv) * Math.PI);
					var y5 = cy + (centre ? radius * 2 / 3 : radius + 8) * Math.sin(-(0.5 - 2 * tv) * Math.PI);
					svg_sm.setAttribute("d", "M " + x4 + " " + y4 + " L " + x5 + " " + y5 + " L 520 " + y5);
					svg_sl.setAttribute("x", 520); 
					svg_sl.setAttribute("y", y5 + 4);
					render_text_svg(svg_sl, "\uEEE8" + (pie[h].value * 100).toFixed(1) + " %" +
							(pie.length === 4 ? " / " + ((pie[h].value / (pie[h].value + pie[(h + 2) % 4].value)) * 100).toFixed(1) + " %  " : "  "));
				}
			}
			
			sigma_value += pie[h].value;
		}
	}
	
	for(t = curves.length - 1; t >= 0; t--)
	{
		ele_id = "svg_curve" + t;
		var svg_curve = document.getElementById(ele_id);
		var svg_curve_points = "";
		
		if(!curves[t] || !curves[t].length)
		{
			//if(svg_curve) remove_from_parent(svg_curve);
			if(svg_curve) svg_curve.style.visibility = "hidden";
			continue;
		}
		
		if(!svg_curve)
		{
			svg_curve = document.createElementNS(svg_xmlns, "polyline");
			svg_curve.id = ele_id;
			if(labels[t].includes("\uEEECK"))
				svg_curve.style.strokeDasharray="6, 6"
			figure_curvegroup.appendChild(svg_curve);
		}
		
		svg_curve.style.visibility = "visible";
		
		// p = transform(curves[t][0]);
		
		for(h = 0; h < curves[t].length; h++)
		{
			p = transform(curves[t][h]);
			
			if(Number.isNaN(p.x) || Number.isNaN(p.y) || !Number.isFinite(p.x) || !Number.isFinite(p.y)) continue;
			
			svg_curve_points = svg_curve_points + (ca_left + p.x).toFixed(2) + "," + (ca_top + p.y).toFixed(2) + " ";
		}
		
		svg_curve.style.stroke = colours[t];
		svg_curve.setAttribute("points", svg_curve_points);
	}
	
	for(t = 0; t < datapoints.length; t++)
	{
		var q = {x: 0, y: 0};
		
		if(xaxistype === axistype_log)
			q.x = Math.log10(datapoints[t].x)
		else
			q.x = datapoints[t].x;
		
		if(yaxistype === axistype_log)
			q.y = Math.log10(datapoints[t].y);
		else
			q.y = datapoints[t].y;
		
		var p = transform(q);
		
		ele_id = "svg_datapoint" + t;
		var svg_cross = document.getElementById(ele_id);
		
		if(!svg_cross)
		{
			svg_cross = document.createElementNS(svg_xmlns, "path");
			svg_cross.id = ele_id;
			svg_cross.style.stroke = "black";
			svg_cross.style.strokeWidth = 1;
			svg_cross.style.fill = "black";
			figure_datapointgroup.appendChild(svg_cross);
		}
		
		var cross_radius = 3;
		
		svg_cross.setAttribute("d", "M " + (ca_left + p.x - cross_radius) + " " + (ca_top + p.y - cross_radius) + " " +
			"Q " + (ca_left + p.x) + " " + (ca_top + p.y) + " " + (ca_left + p.x - cross_radius) + " " + (ca_top + p.y + cross_radius) + " " +
			"Q " + (ca_left + p.x) + " " + (ca_top + p.y) + " " + (ca_left + p.x + cross_radius) + " " + (ca_top + p.y + cross_radius) + " " +
			"Q " + (ca_left + p.x) + " " + (ca_top + p.y) + " " + (ca_left + p.x + cross_radius) + " " + (ca_top + p.y - cross_radius) + " " +
			"Q " + (ca_left + p.x) + " " + (ca_top + p.y) + " " + (ca_left + p.x - cross_radius) + " " + (ca_top + p.y - cross_radius));
		
	}
	
	for(t = 0; t < curves.length; t++)
	{
		var nodraw = false;
		
		if(curves[t] === undefined || labels[t] === undefined || labels[t] === "") nodraw = true;
		
		if(legends[t] || legends[t] === 0)
		{
			ele_id = "svg_legendrect" + t;
			var svg_rect = document.getElementById(ele_id);
			
			ele2_id = "svg_legendlabel" + t;
			var svg_label = document.getElementById(ele2_id);
			
			if(nodraw)
			{
				if(svg_rect) remove_from_parent(svg_rect);
				if(svg_label) remove_from_parent(svg_label);
				continue;
			}
			
			if(!svg_rect)
			{
				svg_rect = document.createElementNS(svg_xmlns, "rect");
				svg_rect.id = ele_id;
				svg_rect.style.stroke = "none";
				figure_legendgroup.appendChild(svg_rect);
			}
			
			svg_rect.style.fill = colours[t];
			svg_rect.setAttribute("x", ca_left + ca_width + 7);
			svg_rect.setAttribute("y", ca_top + legends[t] * 40);
			svg_rect.setAttribute("width", cwidth - ca_left - ca_width - 14);
			svg_rect.setAttribute("height", 30);
			
			if(!svg_label)
			{
				svg_label = document.createElementNS(svg_xmlns, "text");
				svg_label.id = ele2_id;
				svg_label.style.textAnchor = "middle";
				svg_label.style.fill = "white";
				figure_legendgroup.appendChild(svg_label);
			}
			
			render_text_svg(svg_label, labels[t]);
			svg_label.setAttribute("x", ca_left + ca_width + 5 + 0.5 * (cwidth - ca_left - ca_width - 10));
			svg_label.setAttribute("y", ca_top + legends[t] * 40 + 18);
		}
		else
		{
			if(curves[t]) p = transform(curves[t][curves[t].length - 1]);
			
			ele_id = "svg_curvelabel" + t;
			var svg_label = document.getElementById(ele_id);
			
			if(nodraw || p.y > ca_height + 4)
			{
				if(svg_label) remove_from_parent(svg_label);
				continue;
			}
			
			if(!svg_label)
			{
				svg_label = document.createElementNS(svg_xmlns, "text");
				svg_label.id = ele_id;
				figure_legendgroup.appendChild(svg_label);
			}
			
			render_text_svg(svg_label, labels[t]);
			svg_label.style.fill = colours[t];
			svg_label.setAttribute("x", ca_left + Math.min(p.x, ca_width));
			svg_label.setAttribute("y", ca_top + p.y - 2);
		}
	}
	
	if(cwidth - ca_left - ca_width > 0)
	{
		ele_id = "svg_mask1";
		var svg_rect = document.getElementById(ele_id);
		
		if(!svg_rect)
		{
			svg_rect = document.createElementNS(svg_xmlns, "rect");
			svg_rect.id = ele_id;
			svg_rect.style.fill = "white";
			svg_rect.style.stroke = "none";
			figure_overlaygroup.appendChild(svg_rect);
		}
		
		svg_rect.setAttribute("x", ca_left + ca_width + 2);
		svg_rect.setAttribute("y", 0);
		svg_rect.setAttribute("width", cwidth - ca_left - ca_width + 1);
		svg_rect.setAttribute("height", cheight);
	}
	
	if(cheight - ca_top - ca_height > 0)
	{
		ele_id = "svg_mask2";
		var svg_rect = document.getElementById(ele_id);
		
		if(!svg_rect)
		{
			svg_rect = document.createElementNS(svg_xmlns, "rect");
			svg_rect.id = ele_id;
			svg_rect.style.fill = "white";
			svg_rect.style.stroke = "none";
			figure_overlaygroup.appendChild(svg_rect);
		}
		
		svg_rect.setAttribute("x", 0);
		svg_rect.setAttribute("y", ca_top + ca_height);
		svg_rect.setAttribute("width", cwidth);
		svg_rect.setAttribute("height", cheight - ca_top - ca_height + 3);
	}
	
	ele_id = "svg_axisline";
	var svg_axes = document.getElementById(ele_id);
	
	if(!svg_axes)
	{
		svg_axes = document.createElementNS(svg_xmlns, "polyline");
		svg_axes.id = ele_id;
		figure_linegroup.appendChild(svg_axes);
	}
	
	svg_axes.setAttribute("points", ""
		+ ca_left + ",0 "
		+ ca_left + "," + (ca_top + ca_height) + " "
		+ cwidth + "," + (ca_top + ca_height)
	);
	
	var xtickmarks = 10;
	var ytickmarks = 10;
	
	var xmagnitude, xstep;
	var ymagnitude, ystep;
	
	switch(xaxistype)
	{
		case axistype_lin:
		{
			xmagnitude = Math.floor(Math.log10(xmax - xmin));
			xstep = (xmax - xmin) / xtickmarks;
			
			break;
		}
		case axistype_log:
		{
			xstep = (xmax - xmin) / Math.abs(xmax - xmin);
			xtickmarks = Math.round(Math.abs(xmax - xmin));
			
			break;
		}
	}
	
	switch(yaxistype)
	{
		case axistype_lin:
		{
			ymagnitude = Math.floor(Math.log10(ymax - ymin));
			ystep = (ymax - ymin) / ytickmarks;
			
			break;
		}
		case axistype_log:
		{
			ystep = (ymax - ymin) / Math.abs(ymax - ymin);
			ytickmarks = Math.round(Math.abs(ymax - ymin));
			
			break;
		}
	}
	
	var labelstr = "";
	
	var xlist = document.getElementsByClassName("xtick");
	for(h = 0; h < xlist.length;)
		if(Number(xlist[h].id.substr(-2)) > xtickmarks + 10)
			remove_from_parent(xlist[h]);
		else
			h++;
	
	var ylist = document.getElementsByClassName("ytick");
	for(h = 0; h < ylist.length;)
		if(Number(ylist[h].id.substr(-2)) > ytickmarks + 10)
			remove_from_parent(ylist[h]);
		else
			h++;
	
	for(h = 0; h <= Math.max(xtickmarks, ytickmarks); h++)
	{
		var tp = {
			x: xmin + (h / xtickmarks) * (xmax - xmin),
			y: ymin + (h / ytickmarks) * (ymax - ymin)
		};
		
		p = transform(tp);
		
		if(h <= xtickmarks)
		{
			switch(xaxistype)
			{
				case axistype_lin:
				{
					labelstr = replace_minus_signs(((xmin + h * xstep) / Math.pow(10, xmagnitude)).toFixed(2));
					break;
				}
				case axistype_log:
				{
					if(xmin + h * xstep === 0)
						labelstr = "1";
					else
						labelstr = replace_minus_signs("10\uEEE1" + (xmin + h * xstep).toString());
					break;
				}
			}
			
			ele_id = "svg_xtick" + (h + 10);
			var svg_tick = document.getElementById(ele_id);
			
			if(!svg_tick)
			{
				svg_tick = document.createElementNS(svg_xmlns, "line");
				svg_tick.id = ele_id;
				svg_tick.setAttribute("class", "xtick");
				figure_linegroup.appendChild(svg_tick);
			}
			
			svg_tick.setAttribute("x1", ca_left + p.x);
			svg_tick.setAttribute("y1", ca_top + ca_height);
			svg_tick.setAttribute("x2", ca_left + p.x);
			svg_tick.setAttribute("y2", ca_top + ca_height + 4);
			
			ele_id = "svg_xlabel" + (h + 10);
			var svg_label = document.getElementById(ele_id);
			
			if(!svg_label)
			{
				svg_label = document.createElementNS(svg_xmlns, "text");
				svg_label.id = ele_id;
				svg_label.setAttribute("text-anchor", "middle");
				svg_label.setAttribute("class", "xtick");
				figure_overlaygroup.appendChild(svg_label);
			}
			
			render_text_svg(svg_label, labelstr);
			svg_label.setAttribute("x", ca_left + p.x);
			svg_label.setAttribute("y", ca_top + ca_height + 16);
		}
		
		if(h <= ytickmarks)
		{
			switch(yaxistype)
			{
				case axistype_lin:
				{
					labelstr = replace_minus_signs(((ymin + h * ystep) / Math.pow(10, ymagnitude)).toFixed(2));
					break;
				}
				case axistype_log:
				{
					if(ymin + h * ystep === 0)
						labelstr = "1";
					else
						labelstr = replace_minus_signs("10\uEEE1" + (ymin + h * ystep).toString());
					break;
				}
			}
			
			ele_id = "svg_ytick" + (h + 10);
			var svg_tick = document.getElementById(ele_id);
			
			if(!svg_tick)
			{
				svg_tick = document.createElementNS(svg_xmlns, "line");
				svg_tick.id = ele_id;
				svg_tick.setAttribute("class", "ytick");
				figure_linegroup.appendChild(svg_tick);
			}
			
			svg_tick.setAttribute("x1", ca_left);
			svg_tick.setAttribute("y1", ca_top + p.y);
			svg_tick.setAttribute("x2", ca_left - 4);
			svg_tick.setAttribute("y2", ca_top + p.y);
			
			ele_id = "svg_ylabel" + (h + 10);
			var svg_label = document.getElementById(ele_id);
			
			if(!svg_label)
			{
				svg_label = document.createElementNS(svg_xmlns, "text");
				svg_label.id = ele_id;
				svg_label.setAttribute("text-anchor", "end");
				svg_label.setAttribute("class", "ytick");
				figure_overlaygroup.appendChild(svg_label);
			}
			
			render_text_svg(svg_label, labelstr);
			svg_label.setAttribute("x", ca_left - 8);
			svg_label.setAttribute("y", ca_top + p.y + 4);
		}
	}
	
	var xlabelstr = "", ylabelstr = "";
	
	function magnitude_string(magnitude, no_trailing_space)
	{
		if(!magnitude)
			return "";
		else
			return "10\uEEE1" + replace_minus_signs(magnitude.toString()) + "\uEEE0" + (no_trailing_space === true ? "" : " ");
	}
	
	switch(appmode)
	{
		case appmode_ligand:
		{
			xlabelstr = (xscale_alternative ? "Free" : "Total") + " ligand concentration (" + magnitude_string(xmagnitude) + "mol l\uEEE1\u22121\uEEE0)";
			ylabelstr = scale_absolute ? "Protein species concentration (mol l\uEEE1\u22121\uEEE0)" : "Relative protein species amount";
			break;
		}
		case appmode_homodimer:
		{
			xlabelstr = (xscale_alternative ? "Free" : "Total") + " protein concentration (" + magnitude_string(xmagnitude) + "mol l\uEEE1\u22121\uEEE0)";
			ylabelstr = scale_absolute ? "Protein species concentration (mol l\uEEE1\u22121\uEEE0)" : "Relative protein species amount";
			break;
		}
		case appmode_ligands:
		{
			xlabelstr = (xscale_alternative ? "Total ligand L" : "Total protein P") + " concentration (mol l\uEEE1\u22121\uEEE0)";
			ylabelstr = scale_absolute ? "Protein species concentration (mol l\uEEE1\u22121\uEEE0)" : "Relative protein species amount";
			break;
		}
		case appmode_receptors:
		{
			xlabelstr = (xscale_alternative ? "Free" : "Total") + " ligand L concentration (" + magnitude_string(xmagnitude) + "mol l\uEEE1\u22121\uEEE0)";
			// ylabelstr = scale_absolute ? "Protein species concentration (mol l\uEEE1\u22121\uEEE0)" : "Protein species concentration (" + magnitude_string(ymagnitude) + "mol l\uEEE1\u22121\uEEE0)";
			switch(scale_absolute)
			{
				case 2: ylabelstr = "Receptor P specificity" + (ymagnitude ? " / " + magnitude_string(ymagnitude, true) : ""); break;
				case 1: ylabelstr = "Protein species concentration (mol l\uEEE1\u22121\uEEE0)"; break;
				case 0: ylabelstr = "Protein species concentration (" + magnitude_string(ymagnitude) + "mol l\uEEE1\u22121\uEEE0)"; break;
			}
			break;
		}
	}
	
	ele_id = "svg_xaxislabel";
	var svg_label = document.getElementById(ele_id);
	
	if(!svg_label)
	{
		svg_label = document.createElementNS(svg_xmlns, "text");
		svg_label.id = ele_id;
		svg_label.setAttribute("text-anchor", "end");
		svg_label.setAttributeNS(xml_xmlns, "xml:space", "preserve");
		figure_overlaygroup.appendChild(svg_label);
	}
	
	render_text_svg(svg_label, xlabelstr);
	svg_label.setAttribute("x", cwidth - 4);
	svg_label.setAttribute("y", ca_top + ca_height + 36);
	
	ele_id = "svg_yaxislabel";
	var svg_label = document.getElementById(ele_id);
	
	if(!svg_label)
	{
		svg_label = document.createElementNS(svg_xmlns, "text");
		svg_label.id = ele_id;
		svg_label.setAttribute("text-anchor", "end");
		svg_label.setAttributeNS(xml_xmlns, "xml:space", "preserve");
		figure_overlaygroup.appendChild(svg_label);
	}
	
	render_text_svg(svg_label, ylabelstr);
	svg_label.setAttribute("transform", "matrix(0 -1 1 0 " + (ca_left - 44) + " 4)");
	svg_label.setAttribute("x", 0);
	svg_label.setAttribute("y", 0);
}
