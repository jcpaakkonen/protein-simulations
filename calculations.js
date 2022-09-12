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

var datalabels = [];

datalabels[appmode_ligand] = ["[PL]", "[P]"];
datalabels[appmode_homodimer] = ["[P<sub>2</sub>]", "[P]"];
datalabels[appmode_ligands] = ["[PL]", "[PL\u2032]", "[P]"];
datalabels[appmode_receptors] = ["[PL]", "[P\u2032L]", "[P]", "[P\u2032]"];

function calculate_ligand_free(E_0, S, K_D)
{
	var v = E_0 * S / (S + K_D);
	return { d: [v, E_0 - v], m: [1, 1] };
}

function calculate_ligand_total(E_0, S_0, K_D)
{
	/*
	// This the calculation described in the supporting information.
	// It is numerically unstable: when [L] -> 0 and [P] -> E_0,
	// precision is lost and [PL] becomes noisy.
	var b = (E_0 - S_0 + K_D) / 2;
	var L = -b + Math.sqrt(b * b + S_0 * K_D);
	var P = E_0 * K_D / (L + K_D);
	return { d: [E_0 - P, P, L], m: [1, 1, 0] };
	*/
	
	// This is the equivalent calculation but numerically stable.
	// When [PL] is calculated first, precision is never significantly lost
	// with possible parameter values in the simulation.
	var minus_b = (S_0 + E_0 + K_D) / 2;
	var v = minus_b - Math.sqrt(minus_b * minus_b - S_0 * E_0);
	return { d: [v, E_0 - v, S_0 - v], m: [1, 1, 0] };
}

function calculate_homodimer_free(P, K_D)
{
	return { d: [P * P / K_D, P], m: [2, 1] };
}

function calculate_homodimer(E_0, K_D)
{
	var P = (-K_D + Math.sqrt(K_D * K_D + 8 * K_D * E_0)) / 4;
	return { d: [P * P / K_D, P], m: [2, 1] };
}

function calculate_ligands_free(A, c_B, c_C, K_D, K_D2)
{
	var AB = c_B * A / (A + K_D);
	var AC = c_C * A / (A + K_D2);
	var B = c_B - AB;
	var C = c_C - AC;
	
	if(appmode === appmode_ligands)
		return { d: [AB, AC, A, B, C], m: [1, 1, 1, 0, 0] };
	else
		return { d: [AB, AC, B, C, A, AB / AC], m: [1, 1, 1, 1, 0, 0] };
}

function calculate_ligands(c_A, c_B, c_C, K_D, K_D2)
{
	var t1 = 1;
	var t2 = K_D2 + K_D + c_C + c_B - c_A
	var t3 = K_D2 * K_D + c_C * K_D + c_B * K_D2 - c_A * K_D2  - c_A * K_D;
	var t4 = -c_A * K_D2 * K_D;
	
	var A, B, C, AB, AC;
	
	A = solve_cubic_newton(t1, t2, t3, t4, c_A);
	
	if(!(A >= 0 && A < c_A))
		A = solve_cubic_bisection(t1, t2, t3, t4, 0, c_A);
	
	return calculate_ligands_free(A, c_B, c_C, K_D, K_D2);
}

function sum_squared_residuals(arr1, arr2, logarithmic)
{
	var i, r, sum = 0;
	
	if(arr1.length !== arr2.length) return NaN;
	
	if(logarithmic)
	{
		for(i = 0; i < arr1.length; i++)
		{
			// r = Math.log(arr1[i].y) - Math.log(arr2[i].y);
			r = Math.log(arr1[i].y / arr2[i].y);
			sum += r * r;
		}
	}
	else
	{
		for(i = 0; i < arr1.length; i++)
		{
			r = arr1[i].y - arr2[i].y;
			sum += r * r;
		}
	}
	
	return sum;
}

function calculate_lsf()
{
	if(datapoints.length < 1)
	{
		var label = document.getElementById("calcstatus");
		if(label)
		{
			label.innerHTML = "There are no data points.";
			label.style.color = "red";
		}
		
		return;
	}
	
	var k1, k2, k3, k4;
	var m1, m2, m3;
	var k1p, k2p, k3p, k4p;
	var m1p, m2p, m3p;
	
	var k1min, k1max, k1step;
	var k2min, k2max, k2step;
	var k3min, k3max, k3step;
	var k4min, k4max, k4step;
	
	var srsum_min;
	var i, j, n;
	var co;
	
	var calcmode;
	
	var logresid = extmode && document.getElementById("ext_calcmode_log").checked;
	
	for(i = 1; i <= 3; i++)
	{
		var ele = document.getElementById("calcmode" + i);
		
		if(ele && ele.checked)
		{
			calcmode = i;
			break;
		}
	}
	
	var theor = [];
	
	switch(appmode)
	{
		case appmode_ligand:
		case appmode_homodimer:
		case appmode_ligands:
		case appmode_receptors:
		{
			var slider_indices = [];
			var sliders = [];
			var kmin = [];
			var kmax = [];
			var kstep = [];
			var kp = [];
			var k = [];
			var m = [];
			var kp = [];
			var fun;
			
			function fun_ligand_free()
			{
				return calculate_ligand_free(m[5], datapoints[i].x, m[7]);
			}
			
			function fun_ligand_total()
			{
				return calculate_ligand_total(m[5], datapoints[i].x, m[7]);
			}
			
			function fun_homodimer_free()
			{
				return calculate_homodimer_free(datapoints[i].x, m[7]);
			}
			
			function fun_homodimer_total()
			{
				return calculate_homodimer(datapoints[i].x, m[7]);
			}
			
			function fun_ligands()
			{
				return calculate_ligands(datapoints[i].x, m[10], m[3], m[7], m[9]);
			}
			
			function fun_ligands_alt()
			{
				return calculate_ligands(m[5], datapoints[i].x, m[3], m[7], m[9]);
			}
			
			function fun_ligands_free()
			{
				return calculate_ligands_free(datapoints[i].x, m[10], m[3], m[7], m[9]);
			}
			
			switch(appmode)
			{
				case appmode_ligand:
				{
					if(xscale_alternative)
						fun = fun_ligand_free;
					else
						fun = fun_ligand_total;
					break;
				}
				case appmode_homodimer:
				{
					if(xscale_alternative)
						fun = fun_homodimer_free;
					else
						fun = fun_homodimer_total;
					break;
				}
				case appmode_ligands:
				{
					if(xscale_alternative)
						fun = fun_ligands_alt;
					else
						fun = fun_ligands;
					break;
				}
				case appmode_receptors:
				{
					if(xscale_alternative)
						fun = fun_ligands_free;
					else
						fun = fun_ligands;
					break;
				}
				default:
				{
					break;
				}
			}
			
			var fixvallist = document.getElementsByClassName("fixval");
			
			for(i = 0; i < fixvallist.length; i++)
			{
				var index = Number(fixvallist[i].id.substring(6));
				var s = document.getElementById("slider" + index);
				
				m[index] = expval(Number(s.value), expparams[index][0], expparams[index][1], expparams[index][2]);
				
				if(!(fixvallist[i].disabled || !fixvallist[i].checked))
				{
					slider_indices.push(index);
					sliders.push(s);
				}
			}
			
			if(sliders.length < 1)
			{
				var label = document.getElementById("calcstatus");
				if(label)
				{
					label.innerHTML = "There are no free parameters to optimise.";
					label.style.color = "red";
				}
				
				return;
			}
			
			if(sliders.length > datapoints.length)
			{
				var label = document.getElementById("calcstatus");
				if(label)
				{
					label.innerHTML = "The number of data points is smaller than the number of free parameters.";
					label.style.color = "red";
				}
				
				return;
			}
			
			if(datalabels[appmode].length > 0)
			{
				for(co = 0; co < datalabels[appmode].length + (appmode === appmode_receptors ? 1 : 0); co++)
				{
					if(document.getElementById("calcoption" + co).checked) break;
				}
			}
			else
			{
				co = 0;
			}
			
			if(appmode === appmode_receptors && co === 4) co = 5;
			
			for(j = 0; j < (calcmode === 1 ? 2 : 1); j++)
			{
				if(j === 0)
				{
					for(i = 0; i < sliders.length; i++)
					{
						switch(calcmode)
						{
							case 1:
								kmin[i] = Number(sliders[i].min);
								kmax[i] = Number(sliders[i].max);
								kstep[i] = 10;
								break;
							case 2:
								kmin[i] = Number(sliders[i].min);
								kmax[i] = Number(sliders[i].max);
								kstep[i] = 1;
								break;
							case 3:
								kmin[i] = Math.max(Number(sliders[i].min), Number(sliders[i].value) - 5);
								kmax[i] = Math.min(Number(sliders[i].max), Number(sliders[i].value) + 5);
								kstep[i] = 1;
								break;
						}
					}
				}
				else
				{
					// calcmode is always 1 in this case
					
					for(i = 0; i < sliders.length; i++)
					{
						kmin[i] = Math.max(Number(sliders[i].min), kp[i] - 20);
						kmax[i] = Math.min(Number(sliders[i].max), kp[i] + 20);
						kstep[i] = 1;
					}
				}
				
				for(i = 0; i < sliders.length; i++)
				{
					k[i] = kmin[i];
				}
				
				var loop = true;
				
				srsum_min = Number.MAX_VALUE;
				
				while(loop)
				{
					for(i = 0; i < sliders.length; i++)
					{
						m[slider_indices[i]] = expval(k[i],
							expparams[slider_indices[i]][0],
							expparams[slider_indices[i]][1],
							expparams[slider_indices[i]][2]);
					}
					
					for(i = 0; i < datapoints.length; i++)
					{
						var cd = fun();
						var divisor;
						
						if(scale_absolute || appmode === appmode_receptors)
							divisor = 1;
						else
							divisor = cd.d.reduce(function(t,v,i){ return t + v * cd.m[i]; }, 0);
						
						theor[i] = {
							x: datapoints[i].x,
							y: cd.d[co] / divisor * (scale_absolute || appmode === appmode_receptors ? 1 : cd.m[co])
						};
					}
					
					var srsum = sum_squared_residuals(datapoints, theor, logresid);
					
					if(srsum < srsum_min)
					{
						srsum_min = srsum;
						kp = k.slice();
					}
					
					for(i = 0; i < sliders.length; i++)
					{
						if((k[i] += kstep[i]) <= kmax[i])
							break;
						else if(i < sliders.length - 1)
							k[i] = kmin[i];
						else
							loop = false;
					}
				}
			}
			
			var equal_to_initial = true;
			
			if(calcmode === 3)
			{
				for(i = 0; i < sliders.length; i++)
				{
					if(kp[i] !== Number(sliders[i].value))
					{
						equal_to_initial = false;
						break;
					}
				}
			}
			
			for(i = 0; i < sliders.length; i++)
			{
				sliders[i].value = kp[i].toString();
				slider_input(slider_indices[i], true);
			}
			
			if(!equal_to_initial) // only in calcmode 3
			{
				calculate_lsf();
				return;
			}
			
			var label = document.getElementById("calcstatus");
			
			if(label)
			{
				label.innerHTML = "Calculation finished.";
				label.style.color = "green";
			}
			
			update();
			return;
		}
		default:
		{
			break;
		}
	}
}

function solve_cubic_newton(q1, q2, q3, q4, a0)
{
	// Cubic equation is solved using the Newton method.
	// This can fail to converge to the correct solution depending
	// on the choice of the initial guess.
	//
	// Equation: q1 * a^3 + q2 * a^2 + q3 * a + q4 = 0
	// initial guess: a0
	
	var a2, a1 = a0, cycles = 0;
	
	function val(x)
	{
		return (q1*x*x*x + q2*x*x + q3*x + q4);
		// return (q1*Math.pow(x,3) + q2*Math.pow(x,2) + q3*x + q4);
	}
	
	function der(x)
	{
		return (3*q1*x*x + 2*q2*x + q3);
		// return (3*q1*Math.pow(x,2) + 2*q2*x + q3);
	}
	
	// Abort if the derivative is very small because
	// the search would likely not converge
	
	if(val(a1) / der(a1) > 10) return NaN;
	
	while(1)
	{
		a2 = a1 - val(a1) / der(a1);
		
		if(Math.abs(a2 / a1 - 1) < 1e-3) return a2;
		a1 = a2;
		
		// Abort if too many cycles have been run
		// (the search will probably not converge)
		
		if(++cycles >= 20) return NaN; // TODO: try to find ways to reduce this
		
		// BTW: nice videos on the chaotic behaviour of the Newton method:
		// https://www.youtube.com/watch?v=-RdOwhmqP5s
		// https://www.youtube.com/watch?v=LqbZpur38nw
	}
}

function solve_cubic_bisection(q1, q2, q3, q4, amin, amax)
{
	// The root of a cubic equation is searched using the bisection method.
	//
	// Equation: q1 * a^3 + q2 * a^2 + q3 * a + q4 = 0
	// A root between amin and amax is searched.
	// Guaranteed to converge as long as a root exists.
	
	function val(x)
	{
		return (q1*x*x*x + q2*x*x + q3*x + q4);
		// return (q1*Math.pow(x,3) + q2*Math.pow(x,2) + q3*x + q4);
	}
	
	var xmin = amin;
	var xmax = amax;
	var ymin = val(xmin);
	var ymax = val(xmax);
	
	if(ymin === 0) return xmin;
	if(ymax === 0) return xmax;
	
	while(1)
	{
		var xmid = 0.5 * (xmin + xmax);
		var ymid = val(xmid);
		
		// This prevents an eternal loop in some edge cases
		if(xmid === xmin || xmid === xmax) return xmid;
		
		if(Math.abs(ymid) === 0)
		{
			return xmid;
		}
		else if(ymin * ymid > 0)
		{
			xmin = xmid;
			ymin = ymid;
		}
		else
		{
			xmax = xmid;
			ymax = ymid;
		}
		
		if(Math.abs(xmax / xmin - 1) < 1e-3) return xmid;
	}
}
