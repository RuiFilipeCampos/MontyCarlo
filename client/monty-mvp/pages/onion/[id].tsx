import type { NextPage } from "next";
import React from "react";
import * as ch from "@chakra-ui/react";
import { NextRouter, useRouter } from "next/router";

import rd3 from "react-d3-library";

import * as d3 from "d3";

const RD3Component = rd3.Component;

function D3Component({ radius_list }) {
  const [nodeToRender, setNodeToRender] = React.useState(null);
  const ref = React.useRef(null);


  const draw = () => {
    var width: number = ref.current ? ref.current.offsetWidth : 0;
    var node = document.createElement("div");

    var svg = d3
      .select(node)
      .append("svg")
      .attr("width", width)
      .attr("height", 300);

    let region_number: number = 1;
    let previous_radius: number = 0;
    for (let _radius of radius_list) {
      let current_radius: number = (_radius / 100) * width;

      let arc = d3.arc()
        .innerRadius(current_radius - 2)
        .outerRadius(current_radius)
        .startAngle(45 * (Math.PI / 180)) //converting from degs to radians
        .endAngle(3); //just radians

      svg
        .append("path")
        .attr("d", arc);

      svg
        .append("text")
        .text(`Region ${region_number}`)
        .attr("y", 25)
        .attr("x", 30 + (radius_list[region_number - 2] / 100) * width);

      previous_radius = current_radius;
      ++region_number;
    }

    svg
      .append("text")
      .text(`Region ${region_number}`)
      .attr("y", 25)
      .attr("x", 30 + (radius_list[region_number - 2] / 100) * width);

    setNodeToRender(node);
  };

  React.useEffect(() => draw(), [radius_list, nodeToRender]);
  React.useEffect(() => window.addEventListener("resize", draw), []);

  return (
    <>
      <ch.Flex ref={ref} w="full" h="300">
        {nodeToRender ? <RD3Component data={nodeToRender} /> : <></>}
      </ch.Flex>
    </>
  );
}

interface Routes {
  previous: string,
  next: string
}



/** Previous and Next Buttons */
const FormNavigation = ({ previous, next }: Routes): JSX.Element => {
  const router: NextRouter = useRouter();

  return (
    <>
      <ch.HStack w="full" align="left">
        <ch.Spacer />
        <ch.Button onClick={() => router.push(previous)}>
          Previous
        </ch.Button>
        {
          !next ? <></> : <ch.Button colorScheme="facebook" onClick={() => router.push(next)}>
            Next
          </ch.Button>
        }

      </ch.HStack>
    </>
  );
};



const MaterialForms = (props: any): JSX.Element[] => {
  let to_return: JSX.Element[] = [<></>];


  const handleChange = (val, i) => {
    props.materialList[i] = val
    props.setMaterialList(props.materialList)
  }

  to_return.push(
    <>
      <ch.FormControl >
        <ch.FormLabel>
          Material Definition
        </ch.FormLabel>
      </ch.FormControl>
    </>
  )

  for (let i = 0; i < props.n; ++i) {
    to_return.push(
      <>
        <ch.FormControl isRequired>
          <ch.HStack>
            <ch.FormLabel>
              {`Region ${i + 1}:`}
            </ch.FormLabel>
            <ch.RadioGroup onChange={(val) => handleChange(val, i)}>
              <ch.HStack spacing={8} w="full">
                <ch.Radio value="air">Air Dry (sea level)</ch.Radio>
                <ch.Radio value="water">Liquid Water (PTN) </ch.Radio>
                <ch.Radio value="gold">Solid Gold</ch.Radio>
              </ch.HStack>
            </ch.RadioGroup>
          </ch.HStack>
        </ch.FormControl>
      </>
    )


  }

  to_return.push(
    <>
      <ch.FormControl isRequired>
        <ch.HStack>
          <ch.FormLabel>
            {`Region ${props.n + 1}:`}
          </ch.FormLabel>
          <ch.RadioGroup value="gold" isDisabled>
            <ch.HStack spacing={8} w="full" >
              <ch.Radio value="air">Air Dry (sea level)</ch.Radio>
              <ch.Radio value="water">Liquid Water (PTN) </ch.Radio>
              <ch.Radio value="gold">Solid Gold</ch.Radio>
            </ch.HStack>
          </ch.RadioGroup>
        </ch.HStack>
      </ch.FormControl>
    </>
  )

  to_return.push(
    <>
      <ch.FormControl >
        <ch.FormHelperText>
          Low density materials like air offer little resistance to
          most kinds of particle radiation, just look around, that clear
          picture in your brain is because your eyes are capturing low energy
          EM radiation that has travelled in a straight path, unimpeded.
          <b>The last region has Gold set to default in order to protect the server</b>. Low density materials
          extended to infinity may cause particles to travel very large distances.
          Thus causing the simulation to either hault or to return unusual coordinate
          values (by exceeding the maximum value allowed by the C double).


        </ch.FormHelperText>
      </ch.FormControl>
    </>
  )


  to_return.push(
    <>
      <br />
      <ch.Text >
      </ch.Text>
    </>
  )

  return to_return;
};

/** Centered Spinner */
const Loading = () => <>
  <ch.Flex w="100vw" h="100vh">
    <ch.HStack w="full" h="full" align="center" >
      <ch.VStack w="full" h="full" align="center" >
        <ch.Spacer />
        <ch.Spinner size='xl' />
        <ch.Spacer />
      </ch.VStack>
    </ch.HStack>
  </ch.Flex>
</>


const Title = () => <>
  <ch.HStack w="full">
    <ch.Heading>
      The Onion
    </ch.Heading>
    <ch.Spacer />
  </ch.HStack>
</>

const Layout = (props) => {
  return <>  <ch.Flex w="100vw" h="full">
    <ch.Flex shadow="md" w="full" h="full" p="10" m="50">
      <ch.VStack w="full">
        <Title />
        <ch.Spacer />
        <ch.Box w="full" h="full">
          <form onSubmit={props.onSubmit}>
            <ch.VStack w="full" h="full" spacing={10}>
              {props.children}


              <ch.FormControl>
                <ch.Button type="submit" colorScheme="facebook" w="full" >
                  Simulate !
                </ch.Button>
              </ch.FormControl>
              <FormNavigation previous="/onion" next={``} />
            </ch.VStack>

          </form>
        </ch.Box>

      </ch.VStack>
    </ch.Flex>
  </ch.Flex>
  </>
}


const Onion: NextPage = (): JSX.Element => {
  // All states for handling the form information.
  const [radiusList, setRadiusList] = React.useState([0]);
  const [materialList, setMaterialList] = React.useState(['']);
  const [energy, setEnergy] = React.useState(0);
  const [particle, setParticle] = React.useState('');

  // Getting the number of shells.
  const router = useRouter();
  const { id } = router.query;

  React.useEffect(() => {

    if (id == undefined) return;

    let n: number = parseInt(id);
    let radiusList0 = [];
    let materialList0 = [];

    for (let i = 0; i < n; ++i) {
      radiusList0.push(10 * (i + 1));
      materialList0.push('none');
    }

    setRadiusList(radiusList0);
    setMaterialList(materialList0)

  }, [id]);

  if (id == undefined) return <Loading />;
  if (parseInt(id) != radiusList.length) return <Loading />;
  // The code block above will loop until `id` is an number.



  const makeRangeSliders = () => {
    let ret: JSX.Element[] = [<></>];
    for (let i = 0; i < radiusList.length; ++i) {
      ret.push(
        <>
          <ch.RangeSliderThumb
            w="35px"
            index={i}
            children={
              <ch.Text fontSize="10">
                {
                  `${radiusList[i]} cm`
                }
              </ch.Text>
            }
          />
        </>
      );
    }

    return ret;
  };

  const handleSubmit = (e: any) => {
    e.preventDefault();

    let parameters = {
      n: id,
      radiusList: radiusList,
      materialList: [...materialList, 'gold'],
      particle: particle,
      energy: energy
    }


    const u = new URLSearchParams(parameters).toString();
    window.open('/onion/simulate' + '?' + u)
  }


  return (
    <>
      <Layout onSubmit={handleSubmit}>
        <ch.Box w="full">



          {/* The range sliders */}
          <ch.FormControl>
            <ch.FormLabel>
              Geometry Definition
            </ch.FormLabel>
            <ch.RangeSlider
              aria-label={["min", "max"]}
              defaultValue={radiusList}
              minStepsBetweenThumbs={10}
              onChange={ /* Prevents values below 10cm. */
                (value: number[]) => {
                  if (value[0] <= 10) value[0] = 10;
                  setRadiusList(value);
                }
              }
            >
              <ch.RangeSliderTrack>
                <ch.RangeSliderFilledTrack />
              </ch.RangeSliderTrack>
              {makeRangeSliders()}
            </ch.RangeSlider>
            <D3Component radius_list={radiusList} />
          </ch.FormControl>
        </ch.Box>



        {/* To choose the materials */}
        <ch.Box w="full">
          <MaterialForms
            n={radiusList.length}
            materialList={materialList}
            setMaterialList={setMaterialList}
          />
        </ch.Box>




        <ch.VStack spacing={5} w="full">


          {/* To choose the particle type */}
          <ch.FormControl isRequired>
            <ch.FormControl >
              <ch.FormLabel>
                Source Definition
              </ch.FormLabel>
            </ch.FormControl>
            <ch.FormLabel>
              Particle Type
            </ch.FormLabel>

            <ch.Select placeholder="Select Particle" onChange={val => setParticle(val)}>
              <option value="photon">
                Photon
              </option>
              <option value="electron">
                Electron
              </option>
              <option value="positron">
                Positron
              </option>
            </ch.Select>
            <ch.FormHelperText>
              Photons don't interact much with matter (generally result in quicker simulations). Charged particles like Positrons and Electrons interact much more with the medium.
            </ch.FormHelperText>
          </ch.FormControl>

          {/* To choose the particle type */}
          <ch.FormControl isRequired>
            <ch.FormLabel>
              Energy
            </ch.FormLabel>
            <ch.Select placeholder="Select Initial Energy" onChange={val => setEnergy(val)}>
              <option value="1">
                1MeV
              </option>
              <option value="10">
                10MeV
              </option>
              <option value="30">
                30MeV
              </option>
            </ch.Select>
            <ch.FormHelperText>
              The higher the energy, the longer the simulation will take.
            </ch.FormHelperText>
          </ch.FormControl>
        </ch.VStack>



      </Layout>
    </>
  );
};

export default Onion;
