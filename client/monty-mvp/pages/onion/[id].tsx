import type { NextPage } from "next";
import Head from "next/head";
import Image from "next/image";
import styles from "../styles/Home.module.css";
import React from "react";
import * as ch from "@chakra-ui/react";
import { useRouter } from "next/router";

import rd3 from "react-d3-library";

import * as d3 from "d3";

const RD3Component = rd3.Component;

function D3Component({ radius }) {
  const [nodeToRender, setNodeToRender] = React.useState(null);

  const ref = React.useRef(null);

  React.useEffect(() => {
    console.log("width", ref.current ? ref.current.offsetWidth : 0);
  }, [ref.current]);
  let dr;

  const draw = () => {
    var width = ref.current ? ref.current.offsetWidth : 0;
    var node = document.createElement("div");
    var svg = d3
      .select(node)
      .append("svg")
      .attr("width", width)
      .attr("height", 300);

    let region_number = 1;
    let previous_radius = 0;
    for (let radius_ of radius) {
      let current_radius = (radius_ / 100) * width;
      let arc = d3
        .arc()
        .innerRadius(current_radius - 2)
        .outerRadius(current_radius)
        .startAngle(45 * (Math.PI / 180)) //converting from degs to radians
        .endAngle(3); //just radians

      svg.append("path").attr("d", arc);
      dr = current_radius - previous_radius;

      svg
        .append("text")
        .text(`Region ${region_number}`)
        .attr("y", 25)
        .attr("x", 30 + (radius[region_number - 2] / 100) * width);

      previous_radius = current_radius;
      ++region_number;
    }

    svg
      .append("text")
      .text(`Region ${region_number}`)
      .attr("y", 25)
      .attr("x", 30 + (radius[region_number - 2] / 100) * width);

    console.log("last", (radius[region_number - 2] / 100) * width);

    setNodeToRender(node);
  };

  React.useEffect(() => draw(), [radius, nodeToRender]);
  React.useEffect(() => window.addEventListener("resize", draw), []);

  return (
    <ch.Flex ref={ref} w="full" h="300">
      {nodeToRender ? <RD3Component data={nodeToRender} /> : <></>}
    </ch.Flex>
  );
}

const FormNavigation = ({ previous, next }) => {
  const router = useRouter();

  return (
    <ch.HStack w="full" align="left">
      <ch.Spacer />
      <ch.Button onClick={() => router.push(previous)}>
          Previous
      </ch.Button>
      <ch.Button colorScheme="facebook" onClick={() => router.push(next)}>
        Next
      </ch.Button>
    </ch.HStack>
  );
};


const Demo: NextPage = () => {
  const [radiusList, setRadiusList] = React.useState([0]);
  const router = useRouter();
  const { id } = router.query;

  React.useEffect(() => {
    if (id == undefined) return;
    let n = parseInt(id);
    let initial_state = [];
    for (let i = 0; i < n; ++i) {
      initial_state.push(10 * (i + 1));
    }

    setRadiusList(initial_state);

  }, [id]);

  if (id == undefined) return <></>;
  if (parseInt(id) != radiusList.length) return <></>;

  const handleSliderChange = (value) => {
    if (value[0] <= 10) value[0] = 10;
    setRadiusList(value);
  };

  const makeList = () => {
    let to_return = [<></>];
    for (let i = 0; i < radiusList.length; ++i) {
      to_return.push(
        <ch.RangeSliderThumb
          w="35px"
          index={i}
          children={<ch.Text fontSize="10">{`${radiusList[i]} cm`}</ch.Text>}
        />
      );
    }

    return to_return;
  };

  const makeMaterialForm = () => {
    let to_return = [<></>];
    for (let i = 0; i < radiusList.length + 1; ++i) {
      to_return.push(
        <ch.FormControl>
          <ch.FormLabel>
            {`Material Definition of Region ${i + 1}`}
          </ch.FormLabel>
          <ch.Select placeholder="Select Material">
            <option value="sphere">Water</option>
            <option value="onion">Gold</option>
            <option value="sphere-cut">Air</option>
          </ch.Select>
        </ch.FormControl>
      );
    }

    return to_return;
  };

  return (
    <>
      <ch.Flex w="100vw" h="full">
        <ch.Flex shadow="md" w="full" h="full" p="10" m="50">
          <ch.VStack w="full">
            <ch.HStack w="full">
              <ch.Heading>The Onion</ch.Heading>
              <ch.Spacer />
            </ch.HStack>
            <ch.Spacer />
            <ch.FormControl>
              <ch.FormLabel>Geometry Definition</ch.FormLabel>

              <ch.RangeSlider
                aria-label={["min", "max"]}
                defaultValue={radiusList}
                onChange={handleSliderChange}
                minStepsBetweenThumbs={10}
              >
                <ch.RangeSliderTrack>
                  <ch.RangeSliderFilledTrack />
                </ch.RangeSliderTrack>
                {makeList()}
              </ch.RangeSlider>
              <D3Component radius={radiusList} />

              {makeMaterialForm()}
            </ch.FormControl>

            <FormNavigation previous="/demo" next={`/onion/3/simulate`} />
          </ch.VStack>
        </ch.Flex>
      </ch.Flex>
    </>
  );
};

export default Demo;
